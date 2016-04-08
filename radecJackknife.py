import numpy as np

class radecJackknife:
    def __init__(self, lower_boundary):
        self.rotation_matrix = np.array([[1.,0],[0.,1.0]])
        self.lower_boundary = lower_boundary
        return None

    def radec_to_xyz(self, ra, dec):
         phi = np.deg2rad(ra)
         theta = np.pi/2. - np.deg2rad(dec)
         sintheta = np.sin(theta)
         x = sintheta * np.cos(phi)
         y = sintheta * np.sin(phi)
         z = np.cos(theta)
         return x, y, z

    def xyz_to_radec(self, x,y,z):
         theta = np.arccos(z)
         phi = np.arctan(y/x)
         dec = np.rad2deg(np.pi/2. - theta)
         ra = np.rad2deg(phi)
         pdb.set_trace()
         return ra, dec

    def rotate_pts(self, ra, dec):
        center_ra = np.mean(ra)
        center_dec = np.mean(dec)
        uu, vv, ww = self.radec_to_xyz(center_ra, center_dec)
        xx, yy, zz = self.radec_to_xyz(ra, dec)
        #see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
        theta = 0.0#np.pi
        temp = (uu*xx+vv*yy+ww*zz)*(1.-np.cos(theta))
        xprime = uu*temp + xx*np.cos(theta)+ (-ww*yy + vv*zz)*np.sin(theta)
        yprime = vv*temp + yy*np.cos(theta)+ ( ww*xx - uu*zz)*np.sin(theta)
        zprime = ww*temp + zz*np.cos(theta)+ (-vv*xx + uu*yy)*np.sin(theta)
        #transform back into original coordiantes
        ra_prime, dec_prime = self.xyz_to_radec(xprime, yprime, zprime)
        pdb.set_trace()
        return (ra_prime, dec_prime)
    
    def remove_discontiguous(self, ra, dec):
        new_ra = np.copy(ra)
        new_dec = np.copy(dec)
        small_ra = np.where(ra < self.lower_boundary)[0]
        new_ra[small_ra] += 360.
        return new_ra, new_dec

    def combine_regions(self, another):
        self.all_boundaries = self.all_boundaries +  another.all_boundaries
        self.nregions = len(self.all_boundaries)        

    def generate_regions(self, ra, dec, nregions, tol, max_iter):
        ra, dec = self.remove_discontiguous(ra, dec)
        #nregions will be somewhat approximate... is this ok?
        npts = len(ra)
        desired_size = np.floor(npts/nregions)
        
        #assume regions is roughly square to get number of declination regions
        num_dec_bins = int(np.sqrt(nregions))
        dec_bins = np.linspace(np.min(dec), np.max(dec), num_dec_bins+1)
        #Adjust dec bins so that each one contains some integer multiple of desired points such that total number of bins is desired number
        #we'll assume contiguous region for now
        dec_hist, dec_bin_temp = np.histogram(dec, bins = dec_bins)
        num_iter = 0
        while (np.max(np.abs(dec_hist/desired_size - np.round(dec_hist/desired_size))) > tol and num_iter < max_iter):
            num_iter += 1
            bin_sizes = dec_bins[1:] - dec_bins[:-1]
            delta_dec_bins = -0.01*bin_sizes*(dec_hist/desired_size - np.round(dec_hist/desired_size))
            dec_bins[1:] += delta_dec_bins
            dec_hist, dec_bin_temp = np.histogram(dec, bins = dec_bins)
            print "dec_hist/desired_size = ", dec_hist/desired_size
            print "******"
        #Adjust bottom and top dec bins to extend farther
        dec_bin_widths = dec_bins[1:]- dec_bins[:-1]
        mean_dec_bin_width = np.mean(dec_bin_widths)
        dec_bins[0] = dec_bins[0] - mean_dec_bin_width
        dec_bins[-1] = dec_bins[-1] + mean_dec_bin_width
        #number of ra bins for each dec bin
        num_ra_bin_list = np.round(dec_hist/desired_size)
        self.all_boundaries = []
        #Now generate ra bins
        for dec_bin_i in xrange(0,len(dec_bins)-1):
            in_dec_bin = np.where((dec > dec_bins[dec_bin_i]) & (dec < dec_bins[dec_bin_i+1]))[0]
            desired_percentiles = np.linspace(0., 100., num_ra_bin_list[dec_bin_i]+1)
            ra_bins = np.percentile(ra[in_dec_bin],desired_percentiles)
            #Adjust left and right ra bins to extend farther
            ra_bin_widths = ra_bins[1:]- ra_bins[:-1]
            mean_ra_bin_width = np.mean(ra_bin_widths)
            ra_bins[0] = ra_bins[0] - mean_ra_bin_width
            ra_bins[-1] = ra_bins[-1] + mean_ra_bin_width
            for ra_bin_i in xrange(0, len(ra_bins)-1):
                self.all_boundaries.append(((dec_bins[dec_bin_i], dec_bins[dec_bin_i+1]), (ra_bins[ra_bin_i],ra_bins[ra_bin_i+1])))
        #Store regions in object
        self.nregions = len(self.all_boundaries)
        
    def label_pts(self, ra, dec):
        ra, dec = self.remove_discontiguous(ra, dec)
        #ra, dec = rotate_pts(ra, dec)
        labels = np.zeros(len(ra))-1.
        for ji in xrange(0,self.nregions):
            min_dec, max_dec = self.all_boundaries[ji][0]
            min_ra, max_ra = self.all_boundaries[ji][1]
            in_jk_region = np.where((ra >= min_ra) & (ra <= max_ra) & (dec >= min_dec) & (dec <= max_dec))[0]
            labels[in_jk_region] = ji
        if (np.min(labels) < 0):
            print "Points not assigned to region!"
            pdb.set_trace()
        return labels
    

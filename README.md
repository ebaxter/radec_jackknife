# radec_jackknife

Given a random catalog (i.e. a set of random points on the sky that appropriately sample the survey volume) this code will generate a set of (approximately) equal sized regions that can be used for the purposes of jackknifing.  The regions lie along lines of constant RA and DEC. Once the regions have been initialized, new points (from the data catalog, for instance) can be quickly assigned to jackknife regions.  

Example:


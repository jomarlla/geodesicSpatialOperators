/**
 * Java implementation of the intersection and point-to-line solutions 
 * for geodesics on the ellipsoid with high accuracy,
 * designed by Jose Martinez-Llario and Sergio Baselga [1]
 *
 * The geodetic theory of the algorithms from Sergio Baselga
 * and  Jose Martinez-Llario is exposed in [2]
 * The algorithm uses the library geographiclib from Charles Karney to
 * perform the direct and inverse geodetic problems [3]
 *
 * [1] New paper (sending status)
 * [2] Baselga, S. & Martinez-Llario, J.C. Stud Geophys Geod (2018) 62: 353. https://doi.org/10.1007/s11200-017-1020-z
 * [3] Karney C.F.F., 2017. GeographicLib, version 1.47. (http://geographiclib.sourceforge.net)
 *
 * Copyright (c) Jose Martinez-Llario (2018-2019) and licensed
 * under the GPL-3.0-or-later License. 
 *
 *  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,e
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 **********************************************************************/

  
/*  Java Methods: STX_4PIntersection and STX_MinDistanceFast 
 *  
 *  Example of use (coordinates in radians):
 *  GeodesicSpatialOp gg = new GeodesicSpatialOp (6378137, 1 / 298.257223563);
 *  gg.STX_4PIntersection (0.506145, 0.733038, -1.343903, 0.680678, 0, 0.104719, -0.383972, 1.117010);
 *  System.out.println("Radians          LonX: " + (gg.lonX) + "  LatX: " + (gg.latX));
 *  
 *  gg.STX_MinDistanceFast (0.506145, 0.733038, -1.343903, 0.680678, -0.383972, 1.117010);
 *  System.out.println("Radians          LonX: " + (gg.lonX) + "  LatX: " + (gg.latX));
 *  
 *  You can use it from the command line as follows:
 *    Calculate the point intersection between geodesics AB and CD
 *      java -jar GeodesicSpatialOp.jar int loA laA loB laB loC laC loD laD wgs84
 *      java -jar GeodesicSpatialOp.jar int 29 42 -77 39 0 6 -22 64 wgs84
 *      
 *    Calculate the minimum distance from geodesic AB to point P
 *      java -jar GeodesicSpatialOp.jar min loA laA loB laB loP laP wgs84
 *      java -jar GeodesicSpatialOp.jar min  29 42 -77 39 -22 64 wgs84
 *      
 *    Run java -jar GeodesicSpatialOp.jar help for more information
 */
 
package jomarlla.geodesy;

import java.util.Date;

import net.sf.geographiclib.Geodesic;
import net.sf.geographiclib.GeodesicData;


public class GeodesicSpatialOp {
	
	private void setAF (double a, double f) {
		this.a = a;
		this.f = f;
		ELLIP = new Geodesic (a, f);		
	}
	
	public GeodesicSpatialOp(double a, double f) {
		setAF (a, f);
	}
	
	public GeodesicSpatialOp() {
		//WGS84 by default
		setAF (6378137, 1 / 298.257223563);
	}
	private double a;
	private double f;
	
	public static double DEG2RAD = Math.PI / 180.0;
	public static double RAD2DEG = 180.0 / Math.PI;
	private boolean sphere;
	private int nIterations = 1;

	private Geodesic ELLIP;
		    
	private double NA;
	private double NB;
	private double MA;
	private double MB;
	
	public double latX;
	public double lonX;

	// Not used function
	// lon and lat in rads
    public void getRadiusNM (double lonA, double latA, double lonB, double latB) {
    	double exc2 = 2*f-f*f;
    	
    	double slatA = Math.sin(latA); 
    	double slatA2 = slatA * slatA;
   	    double slatB = Math.sin(latB); 
   	    double slatB2 = slatB * slatB;
   	    
   	    NA = a/Math.sqrt(1-exc2*slatA2); // % radius of curvature in the prime vertical
        NB = a/Math.sqrt(1-exc2*slatB2); 
    	
        MA = a*(1-exc2)/ Math.pow ( (1-exc2*slatA2), 1.5); // % radius of the meridian ellipse
        MB = a*(1-exc2)/ Math.pow ( (1-exc2*slatB2), 1.5);
	}
    

    public static double normAzRad (double az) {
    	// returns a [-180, 180] range
    	double pi = Math.PI;
    	double pi2 = 2 * Math.PI;
    	
    	while (az > pi2) az -= pi2;
    	while (az < -pi2) az += pi2;
    	
    	if (az < - pi) az += pi2;
    	if (az > + pi) az -= pi2;
    	
    	return az;
    }
    
    public void changeToSphere () {
    	//We never make the calculatios using a local sphere (it is just sperical trigonometry). It is a test.
		f = 0;							
		ELLIP = new Geodesic (a, f);   	    	
    }
    
    //The parameters lon and lat are in rads
    //returns: 
    //  the number of iterations needed
    //  the public members lonX and latX contain the intersection point
    public int
    STX_4PIntersection (double lonA, double latA, double lonB, double latB,
     		            double lonC, double latC,   double lonD, double latD) {
    	   
    	//We never make the calculatios using a local sphere (it is just sperical trigonometry). It is a test.
		if (sphere) changeToSphere ();
		
    	double threshold = 0.0001; // 0.1 mm
    	
    	double A, C, R, sAX, sCX, angle1, angle2;
    	
    	int iterNo = 0;
    	int iterates = 1;
    	
    	R = a;
    	
    	double sAC, aAC, aCA;
    	double lonA2, latA2, lonC2, latC2;

    	    	
    	while (iterates  == 1) {  
    		LineData AB = STX_Inverse (lonA, latA, lonB, latB);    		
    		LineData CD = STX_Inverse (lonC, latC, lonD, latD);
    		LineData AC = STX_Inverse (lonA, latA, lonC, latC);
    		
    		A = AC.azi1 - AB.azi1;
    		C = CD.azi1 - AC.azi2 + Math.PI;
    		// Azimuth CA = CA.azi1 = AC.azi2 + 180
    		// Azimuth AC = AC.azi1 = CA.azi2 + 180
    		
    		sAC = AC.s12;
    		aAC = AC.azi1;
    		
    		angle1 = Math.atan(Math.sin(sAC/R)/(Math.cos(aAC/R)*Math.cos(A)+(1/Math.tan(C))*Math.sin(A)));
  			sAX = R * angle1; 
  			
  			angle2 = Math.atan(Math.sin(sAC/R)/(Math.cos(aAC/R)*Math.cos(C)+(1/Math.tan(A))*Math.sin(C)));
   			sCX = R * angle2;
   			
   			
   			/*
   			// For obtaining the near-antipodal intersection extending the geodesics uncomment this
   			if (iterNo == 0) {
  				sAX = R * ( angle1 + Math.PI );
   				sCX = R * ( angle2 + Math.PI );    			
   			}
   			 */
    			
    		iterNo++;
    		    		
    		LineData AA2 = STX_Direct (lonA, latA, AB.azi1, sAX);
    		latA2 = AA2.lat2;
    		lonA2 = AA2.lon2;
    		
    		//Comprobar si es necesario, no está en el paper
    		/*
    		if (latC == Math.PI/2 || latC == -Math.PI/2 ) {
    			LineData CA2 = STX_Inverse (lonC, latC, AA2.lon2, AA2.lat2);
    			sCX = CA2.s12;
    		}
    		*/
    		    		
    		LineData CC2 = STX_Direct (lonC, latC, CD.azi1, sCX);
    		latC2 = CC2.lat2;
    		lonC2 = CC2.lon2;
    		
    		/*
    		System.out.println("----- Iteration ---- : " + iterNo);    		
    		System.out.println("sAX: " + sAX + "    sCX: " + sCX);
    		System.out.println("latA2-latC2: " +  (latA2-latC2) + "      lonA2-lonC2: " +  (lonA2-lonC2) );
        	System.out.println("DDMMS LonA2: " + ddmmss(lonA2 * RAD2DEG) + "  LatA2: " + ddmmss(latA2 * RAD2DEG));
        	System.out.println("DDMMS LonC2: " + ddmmss(lonC2 * RAD2DEG) + "  LatC2: " + ddmmss(latC2 * RAD2DEG));
    		*/
    		
    		if ( ( Math.abs(sAX) < threshold && Math.abs(sCX) < threshold) ||
       		     ( Math.abs(latA2-latC2) < 1E-14 && Math.abs(lonA2-lonC2) < 1E-14 ) ) {
    		
    			iterates = 0;    	
				latX = latA2;
				lonX = lonA2;    			
    		} else {
    			latA=latA2;
    			lonA=lonA2;
    		    latC=latC2;
    		    lonC=lonC2;
    		}    	    		    	
    	}
    	
    	return iterNo;    	
    }
    
    private  void setSpheroid (String args[], int nArgs)  {
    	nIterations = 1;
    	sphere = false;
    	
    	// Get the spheroid    	
    	if (args[nArgs].equalsIgnoreCase("wgs84")) {
    		a = 6378137;
    		f = 1 / 298.257223563;
    	} else if (args[nArgs].equalsIgnoreCase("grs80")){
    		a = 6378137;
    		f = 1 / 298.257222101;    		
    	}  else if (args[nArgs].equalsIgnoreCase("sphere")){
    		sphere = true;
    		//If the method does not support local sphere then use wgs84
    		a = 6378137;
    		f = 1 / 298.257223563;
    	} else if (args[nArgs].equalsIgnoreCase("inter")){
    		a = 6378388;
    		f = 1 / 297;    		
    	} else if (args[nArgs].equalsIgnoreCase("custom")){ // Should be the value for axes
    	   	a = Double.parseDouble(args[nArgs + 1]);
    	   	f = Double.parseDouble(args[nArgs + 2]);
    	   	nArgs +=2;
    	}
    	
    	nArgs ++;
    	
    	// Iterations
    	if (args.length > nArgs) {
    		nIterations = Integer.parseInt(args[nArgs]);
    	}    	
    	
    	setAF (a, f);
    	
    }
    
  
    private static void mainMin (String[] args)  {
    	String tipo = "min";
    	
    	if (args.length > 0) tipo = args[0];
    	
       	if (args.length < 8 || args.length > 10) {
    			System.out.println("java -jar GeodesicSpatialOp.jar " + tipo + " loA laA loB laB loP laP custom axis_a eccentricity [nIterations=1]");
    			System.out.println("java -jar GeodesicSpatialOp.jar " + tipo + " loA laA loB laB loP laP esferoid=wgs84,grs80,sphere [nIterations=1]");
    			System.out.println("");
    			System.out.println("Examples (coordinates in decimal degrees):");
    			System.out.println("    java -jar GeodesicSpatialOp.jar " + tipo + " 29 42 -77 39 -22 64 wgs84");
    			System.out.println("    java -jar GeodesicSpatialOp.jar " + tipo + " 29 42 -77 39 -22 64 custom 6378137 298.257223563");
    			System.out.println("    java -jar GeodesicSpatialOp.jar " + tipo + " 29 42 -77 39 -22 64 grs80 1000");
    			System.out.println("       Calculate 1000 iterations and shows the running time");
    			System.out.println("");
    			System.out.println("Running with demo: java -jar GeodesicSpatialOp.jar " + tipo + " 29 42 -77 39 -22 64 wgs84 1000");
    			System.out.println("");
    			
    			args = new String [9];    			
    			args[0] = tipo;
     			args[7] = "wgs84"; args[8] = "1000";
     			
     			//args[7] = "sphere"; args[8] = "1";
     		    			
    			//Case 3    			
    			args[1] = "29"; args[2] = "42"; args[3] = "-70"; args[4] = "-35";
    			args[5] = "-22"; args[6] = "64"; 
    			
    			
    			//Case 2		
    			args[1] = "29"; args[2] = "42"; args[3] = "-77"; args[4] = "39";
    			args[5] = "-22"; args[6] = "64"; 
    			  			
    			
    			// Case 1
    			args[1] = "5"; args[2] = "52"; args[3] = "6"; args[4] = "51.4";
    			args[5] = "5.5"; args[6] = "52"; 
    			/*
    		
    		
    			//Case 4
    			args[1] = "56.3"; args[2] = "76.4"; args[3] = "-129.4"; args[4] = "75.3";
    			args[5] = "-36.7"; args[6] = "36"; 
    			
    			
    			//Case 5
    			args[1] = "50"; args[2] = "70"; args[3] = "80"; args[4] = "70";
    			args[5] = "3"; args[6] = "90"; 
    			
    			
    			//Case 6a (two minimum distances)
    			args[1] = "50"; args[2] = "0"; args[3] = "80"; args[4] = "0";
    			args[5] = "3"; args[6] = "90"; 
    			
    			
    			//Case 6
    			args[1] = "50"; args[2] = "70"; args[3] = "-130"; args[4] = "70";
    			args[5] = "140"; args[6] = "30";
    			
    			
    			//Case 8
    			args[1] = "-42"; args[2] = "40"; args[3] = "63"; args[4] = "65.5";
    			args[5] = "15.438098599510266"; args[6] = "68.64690953749886"; //Intersection point from case 9
    			
    			*/
    			
    		}
        	
        	double xy[] = new double[6];   
        	GeodesicSpatialOp jc = null;
        	
        	try {
    	    	//Get the coordinates
    	    	for (int i = 0; i < 6; i++) {
    	    		xy[i] = Double.parseDouble(args[i+1]) * DEG2RAD;
    	    	}
    	    	
    	    	jc = new GeodesicSpatialOp ();
    	    	jc.setSpheroid(args, 7);
    	    	
    		} catch (Exception e) {
        		System.out.println("There are some erros with the parameters : Try java -jar GeodesicSpatialOp.jar help");
        		System.out.println("--> " + e.getLocalizedMessage());
        		System.exit(1);
        	}
        	
       
        	
        	System.out.println("Using  a spheroid (" + jc.a + "," + jc.f + ")");
        	long t0 = new Date().getTime();
        	int iterNo = 0;
        	        	
        	//the new STX_MinDistance gives a high accuracy
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		iterNo = jc.STX_MinDistance (xy[0],xy[1],xy[2],xy[3],xy[4],xy[5]);
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("The algorithm converges in " + iterNo + " iterations");
        	System.out.println("DDMMSS           LonX: " + ddmmss(jc.lonX * RAD2DEG) + "  LatX: " + ddmmss(jc.latX * RAD2DEG));
        	System.out.println("Radians          LonX: " + (jc.lonX) + "  LatX: " + (jc.latX));
        	System.out.println("Decimal Degrees  LonX: " + (jc.lonX * RAD2DEG) + "  LatX: " + (jc.latX * RAD2DEG));
        	
        	// Test A and B from paper
        	double lonX1 = jc.lonX;
        	double latX1 = jc.latX;
        	LineData AB = jc.STX_Inverse (xy[0],xy[1],xy[2],xy[3]);
        	LineData AX = jc.STX_Inverse (xy[0],xy[1],lonX1, latX1);
        	LineData XP = jc.STX_Inverse (lonX1, latX1,xy[4],xy[5]);
        	System.out.println("Minimum distance: : " + XP.s12); 
        	
        	LineData AX2 = jc.STX_Direct(xy[0],xy[1], AB.azi1, AX.s12);
        	LineData X1X2 = jc.STX_Inverse(lonX1, latX1, AX2.lon2, AX2.lat2);
        	
        	LineData X1P2a = jc.STX_Direct(lonX1, latX1, AX.azi2 + Math.PI / 2, XP.s12);
        	LineData X1P2b = jc.STX_Direct(lonX1, latX1, AX.azi2 - Math.PI / 2, XP.s12);
        	LineData PP2a = jc.STX_Inverse(xy[4],xy[5], X1P2a.lon2, X1P2a.lat2);
        	LineData PP2b = jc.STX_Inverse(xy[4],xy[5], X1P2b.lon2, X1P2b.lat2);
        	double d_PP2 = Math.min(PP2a.s12, PP2b.s12);
        	
        	System.out.println("");
        	System.out.println("Checking Test B");
        	System.out.println("  Distance X1X2: " + X1X2.s12);
        	//System.out.println("  Distance PP2a: " + PP2a.s12);
        	//System.out.println("  Distance PP2b: " + PP2b.s12);
        	System.out.println("  Distance PP2: " + d_PP2);
        	System.out.println("  Test B (meters): Maximum (X1X2, PP2): " + Math.max(X1X2.s12, d_PP2));
        	
        	System.out.println("Checking Test A");
        	double az1 = AB.azi1 - AX.azi1;
        	double az2 = normAzRad (Math.PI / 2 - Math.abs(AX.azi2 + Math.PI - XP.azi1) );
        	
        	//XB.azi1 = XA.azi2, X coult be not on the AB
        	LineData XB = jc.STX_Inverse (lonX1, latX1, xy[2],xy[3]);
        	double az3 = normAzRad (Math.PI / 2 - Math.abs(XB.azi1 - XP.azi1) );
        	
        	System.out.println("  Azimuth aAB-aAX: " + az1);
        	System.out.println("  Azimuth 90 - |aXA - aXP|: " + az2);
        	System.out.println("  Azimuth 90 - |aXB - aXP|: " + az3);
        	System.out.println("  Test A (arc seconds) - Maximum (az1,az2,az3): " + 
        	  Math.max( Math.max(Math.abs(az1), Math.abs(az2)), az3) * RAD2DEG * 3600 );
        	// End of Test A and B from paper
           	 
        
    }
 
    private static void mainDirect (String[] args)  {
       	if (args.length < 6 || args.length > 8) {
    			System.out.println("java -jar GeodesicSpatialOp.jar direct loA laA azAB sAB custom axis_a eccentricity [nIterations=1]");
    			System.out.println("java -jar GeodesicSpatialOp.jar direct loA laA azAB sAB esferoid=wgs84,grs80,localsphere [nIterations=1]");
    			System.out.println("");
    			System.out.println("Examples (coordinates in decimal degrees):");
    			System.out.println("    java -jar GeodesicSpatialOp.jar direct 29 42 30 5000 wgs84");
    			System.out.println("    java -jar GeodesicSpatialOp.jar direct 29 42 30 5000 custom 6378137 298.257223563");
    			System.out.println("    java -jar GeodesicSpatialOp.jar direct 29 42 30 5000 grs80 1000");
    			System.out.println("       Calculate 1000 iterations and shows the running time");
    			System.out.println("");
    			System.out.println("Running with demo: java -jar GeodesicSpatialOp.jar direct 29 42 30 5000 wgs84 1000");
    			System.out.println("");
    			
    			args = new String [7];
    			args[0] = "direct";
    			args[1] = "29"; args[2] = "42"; args[3] = "30"; args[4] = "5000";
    			args[5] = "wgs84"; args[6] = "1000";
    		}
        	
        	double xy[] = new double[4];   
        	GeodesicSpatialOp jc = null;
        	
        	try {
    	    	//Get the coordinates
    	    	for (int i = 0; i < 4; i++) {
    	    		xy[i] = Double.parseDouble(args[i+1]);
    	    		if (i != 3) xy[i] *= DEG2RAD;
    	    	}
    	    	 	    	
    	    	jc = new GeodesicSpatialOp ();
    	    	jc.setSpheroid(args, 5);
    	    	
    		} catch (Exception e) {
        		System.out.println("There are some erros with the parameters : Try java -jar GeodesicSpatialOp.jar help");
        		System.out.println("--> " + e.getLocalizedMessage());
        		System.exit(1);
        	}
        	
       
        	
        	System.out.println("Using  a spheroid (" + jc.a + "," + jc.f + ")");
        	long t0 = new Date().getTime();
        	
        	LineData AB = null;
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		AB = jc.STX_Direct (xy[0],xy[1],xy[2],xy[3]);
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	
        	double aBA = normAzRad (AB.azi2);
        	
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("DDMMSS           LonX: " + ddmmss(jc.lonX * RAD2DEG) + "  LatX: " + ddmmss(jc.latX * RAD2DEG) +  " azBA: " + ddmmss(aBA * RAD2DEG));
        	System.out.println("Radians          LonX: " + (jc.lonX) + "  LatX: " + (jc.latX) + "  azBA: " + (aBA));
        	System.out.println("Decimal Degrees  LonX: " + (jc.lonX * RAD2DEG) + "  LatX: " + (jc.latX * RAD2DEG) + "  azBA: " + (aBA * RAD2DEG));
        	
    	
    }
    

    private static void mainInverse (String[] args)  {
       	if (args.length < 6 || args.length > 8) {
    			System.out.println("java -jar GeodesicSpatialOp.jar inverse loA laA loB laB custom axis_a eccentricity [nIterations=1]");
    			System.out.println("java -jar GeodesicSpatialOp.jar inverse loA laA loB laB esferoid=wgs84,grs80,localsphere [nIterations=1]");
    			System.out.println("");
    			System.out.println("Examples (coordinates in decimal degrees):");
    			System.out.println("    java -jar GeodesicSpatialOp.jar inverse 29 42 -77 39 wgs84");
    			System.out.println("    java -jar GeodesicSpatialOp.jar inverse 29 42 -77 39 custom 6378137 298.257223563");
    			System.out.println("    java -jar GeodesicSpatialOp.jar inverse 29 42 -77 39 grs80 1000");
    			System.out.println("       Calculate 1000 iterations and shows the running time");
    			System.out.println("");
    			System.out.println("Running with demo: java -jar GeodesicSpatialOp.jar inverse 29 42 -77 39 wgs84 1000");
    			System.out.println("");
    			
    			args = new String [7];
    			args[0] = "invserse";
    			args[1] = "29"; args[2] = "42"; args[3] = "-77"; args[4] = "39";
    			args[5] = "wgs84"; args[6] = "1000";
    		}
        	
        	double xy[] = new double[4];   
        	GeodesicSpatialOp jc = null;
        	
        	try {
    	    	//Get the coordinates
    	    	for (int i = 0; i < 4; i++) {
    	    		xy[i] = Double.parseDouble(args[i+1]) * DEG2RAD;
    	    	}
    	    	 	    	
    	    	jc = new GeodesicSpatialOp ();
    	    	jc.setSpheroid(args, 5);
    	    	
    		} catch (Exception e) {
        		System.out.println("There are some erros with the parameters : Try java -jar GeodesicSpatialOp.jar help");
        		System.out.println("--> " + e.getLocalizedMessage());
        		System.exit(1);
        	}
        	
       
        	
        	System.out.println("Using  a spheroid (" + jc.a + "," + jc.f + ")");
        	long t0 = new Date().getTime();
        	
        	LineData AB = null;
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		AB = jc.STX_Inverse(xy[0],xy[1],xy[2],xy[3]);
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	
        	double aAB = normAzRad (AB.azi1);
        	double aBA = normAzRad(AB.azi2);
        	double sAB = AB.s12;
        	
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("Distance AB: " + sAB);
        	System.out.println("DDMMSS           azAB: " + ddmmss(aAB * RAD2DEG) +  " azBA: " + ddmmss(aBA * RAD2DEG));
        	System.out.println("Radians          azAB: " + (aAB) + "  azBA: " + (aBA));
        	System.out.println("Decimal Degrees  azAB: " + (aAB * RAD2DEG) + "  azBA: " + (aBA * RAD2DEG));
        	
    	
    }
     
    public static void main(String[] args)  {

    	if (args.length > 0) {
    		if (args[0].equalsIgnoreCase("int")) { mainInt (args); return; }    		  
    		  else if (args[0].equalsIgnoreCase("min")) { mainMin (args); return; }     		
    		  else if (args[0].equalsIgnoreCase("direct")) { mainDirect (args); return; } 
    		  else if (args[0].equalsIgnoreCase("inverse")) { mainInverse (args); return; } 
    	} 

		System.out.println("Examples (lon and lat in degrees):");
		System.out.println("java -jar GeodesicSpatialOp.jar int loA laA loB laB loC laC loD laD wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar min loA laA loB laB loP laP wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar inverse loA laA loB laB wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar direct loA laA azAB sAB wgs84 [nIterations=1]");
		System.out.println("");
		
		//This line is just for debugging
		//mainMin(args);
    }
    
    private static void mainInt (String[] args)  {	
    	
    	if (args.length < 10 || args.length > 12) {
			System.out.println("java -jar GeodesicSpatialOp.jar int loA laA loB laB loC laC loD laD custom axis_a eccentricity [nIterations=1]");
			System.out.println("java -jar GeodesicSpatialOp.jar int loA laA loB laB loC laC loD laD esferoid=wgs84,grs80,localsphere [nIterations=1]");
			System.out.println("");
			System.out.println("Examples (coordinates in decimal degrees):");
			System.out.println("    java -jar GeodesicSpatialOp.jar int 29 42 -77 39 0 6 -22 64 wgs84");
			System.out.println("    java -jar GeodesicSpatialOp.jar int 29 42 -77 39 0 6 -22 64 custom 6378137 298.257223563");
			System.out.println("    java -jar GeodesicSpatialOp.jar int 29 42 -77 39 0 6 -22 64 grs80 1000");
			System.out.println("       Calculate 1000 iterations and shows the running time");
			System.out.println("");
			System.out.println("Running with demo: java -jar GeodesicSpatialOp.jar int 29 42 -77 39 0 6 -22 64 wgs84 1000");
			System.out.println("");
			
			args = new String [11];
			args[0] = "int";
			args[9] = "wgs84"; args[10] = "1000";
			
			//args[9] = "localsphere"; args[10] = "1";
			
			//Case 1 from paper
			args[1] = "14.5"; args[2] = "54"; args[3] = "14.6"; args[4] = "54.2";
			args[5] = "14.4"; args[6] = "54.1"; args[7] = "14.7"; args[8] = "54";
			
			
			
			//Case 2 from paper
			args[1] = "5"; args[2] = "52"; args[3] = "6"; args[4] = "51.4";
			args[5] = "4.5"; args[6] = "51.5"; args[7] = "5.5"; args[8] = "52";
			
			
			//Case 3 from paper
			args[1] = "29"; args[2] = "42"; args[3] = "-77"; args[4] = "39";
			args[5] = "0"; args[6] = "6"; args[7] = "-22"; args[8] = "64";
			
			
			//Case 4 from paper
			args[1] = "-92"; args[2] = "35"; args[3] = "52"; args[4] = "40";
			args[5] = "20"; args[6] = "-8"; args[7] = "-95"; args[8] = "49";
			/*
			
			//Case 5 from paper
			args[1] = "5"; args[2] = "90"; args[3] = "5"; args[4] = "0";
			args[5] = "-30"; args[6] = "70"; args[7] = "40"; args[8] = "70";
			
			//Case 6 from paper
			args[1] = "-175"; args[2] = "80"; args[3] = "5"; args[4] = "0";
			args[5] = "-30"; args[6] = "60"; args[7] = "40"; args[8] = "80";
			
			//Case 9 from paper
			args[1] = "-42"; args[2] = "40"; args[3] = "63"; args[4] = "65.5";
			args[5] = "-41.8"; args[6] = "40"; args[7] = "62.9"; args[8] = "65.6";
			
			*/		
		}
    	
    	double xy[] = new double[8];
 
    	GeodesicSpatialOp jc = null;
    	
    	try {
	    	//Get the coordinates
	    	for (int i = 0; i < 8; i++) {
	    		xy[i] = Double.parseDouble(args[i+1]) * DEG2RAD;
	    	}
	    	
	    	jc = new GeodesicSpatialOp ();
	    	jc.setSpheroid(args, 9);
	    	
		} catch (Exception e) {
    		System.out.println("There are some erros with the parameters : Try java -jar GeodesicSpatialOp.jar help");
    		System.out.println("--> " + e.getLocalizedMessage());
    		System.exit(1);
    	}
    	

  
    	
    	System.out.println("Using " + (jc.sphere ? "a local spehere" : "a spheroid (" + jc.a + "," + jc.f + ")"));
    	long t0 = new Date().getTime();
    	int iterNo = 0;
    	
    	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
    		iterNo = jc.STX_4PIntersection (xy[0],xy[1],xy[2],xy[3],xy[4],xy[5],xy[6],xy[7]);
    	}
    	    	    	
    	double runningTime = new Date().getTime() - t0;
    	System.out.println("Total Running time in mili seconds: " + runningTime);
    	System.out.println("Number of iterations:" + jc.nIterations);
    	System.out.println("");
    	System.out.println("The algorithm converges in " + iterNo + " iterations");
    	System.out.println("DDMMSS           LonX: " + ddmmss(jc.lonX * RAD2DEG) + "  LatX: " + ddmmss(jc.latX * RAD2DEG));
    	System.out.println("Radians          LonX: " + (jc.lonX) + "  LatX: " + (jc.latX));
    	System.out.println("Decimal Degrees  LonX: " + (jc.lonX * RAD2DEG) + "  LatX: " + (jc.latX * RAD2DEG));
    	
    	
    	//Test A and B
    	double lonX1 = jc.lonX;
    	double latX1 = jc.latX;
    	LineData AB = jc.STX_Inverse (xy[0],xy[1],xy[2],xy[3]);
    	LineData AX = jc.STX_Inverse (xy[0],xy[1],lonX1, latX1);
    	LineData AX2 = jc.STX_Direct(xy[0],xy[1], AB.azi1, AX.s12);
    	LineData AB_X1X2 = jc.STX_Inverse (lonX1, latX1, AX2.lon2, AX2.lat2);
    	
    	LineData CD = jc.STX_Inverse (xy[4],xy[5],xy[6],xy[7]);
    	LineData CX = jc.STX_Inverse (xy[4],xy[5],lonX1, latX1);
    	LineData CX2 = jc.STX_Direct(xy[4],xy[5], CD.azi1, CX.s12);
    	LineData CD_X1X2 = jc.STX_Inverse (lonX1, latX1, CX2.lon2, CX2.lat2);
    	
    	System.out.println("");
    	System.out.println("Checking the geodesic AB");
    	System.out.println("  Azimuth AB: " + AB.azi1);
    	System.out.println("  Azimuth AX: " + AX.azi1);
    	System.out.println("  Azimuth AB-AX: " + (AB.azi1 - AX.azi1));
    	System.out.println("  DistanceAB X1X2: " + AB_X1X2.s12);
    	
    	System.out.println("Checking the geodesic CD");
    	System.out.println("  Azimuth CD: " + CD.azi1);
    	System.out.println("  Azimuth CX: " + CX.azi1);
    	System.out.println("  Azimuth CD-CX: " + (CD.azi1 - CX.azi1));
    	System.out.println("  DistanceCD X1X2: " + CD_X1X2.s12);
    	
    	System.out.println("Test A (arc second)");
    	System.out.println("  Azimuth: " + Math.max(Math.abs(AB.azi1 - AX.azi1), Math.abs(CD.azi1 - CX.azi1)) * RAD2DEG * 3600);
    	System.out.println("Test B (meters)");
    	System.out.println("  Distance: " + Math.max(AB_X1X2.s12, CD_X1X2.s12)); 
    	//End of Test A and B
    }
    
    public static String ddmmss (double g) {
  
    	String ddmmss;
    	int d, m;
    	double t1, s;
    	
    	d = (int) g;  // Truncate the decimals
    	t1 = (g - d) * 60;
    	m = (int) t1;
    	s = (t1 - m) * 60;
    	
    	ddmmss = "" + d + " " + m + " " + String.format("%.9f", s);
    	
    	return ddmmss;
    	
    }
    
    public LineData
    STX_Inverse (double lonA, double latA, double lonB, double latB) {    	
		LineData r = new LineData ();
    	r.lat1 = latA;
    	r.lon1 = lonA;
    	r.lat2 = latB;
    	r.lon2 = lonB;
    	
    	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
    	//Geodesic.Inverse returns azimuths (azi1, azi2) in degrees [-180,180]
		
		r.azi1 = g.azi1 * DEG2RAD; //Azimuth AB on point A
		r.azi2 = g.azi2 * DEG2RAD; //Azimuth AB on point B. Azimuth BA on point B = g.azi2 + 180
		r.s12 = g.s12;   
		
		return r;
    }
    
    public LineData
    STX_Direct (double lonA, double latA, double azAB, double sAB) {  
    	LineData r = new LineData ();
    	r.lat1 = latA;
    	r.lon1 = lonA;
    	r.azi1 = azAB;
    	r.s12 = sAB;
    	
    	GeodesicData g = ELLIP.Direct (latA * RAD2DEG, lonA * RAD2DEG, azAB * RAD2DEG, sAB);    		
		
		r.lat2 = g.lat2 * DEG2RAD;
		r.lon2 = g.lon2 * DEG2RAD;   
		r.azi2 =  g.azi2 * DEG2RAD;
		
		return r;
    }
    
    //The parameters lon and lat are in rads
	//returns: 
	//  the number of iterations needed
	//  the public members lonX and latX contain the closest point on the geodesic AB
	public int
	STX_MinDistance (double lonA, double latA, double lonB, double latB,
	 		            double lonP, double latP) {
		//We do NOT make any calculation using a local sphere (spherical trigonometry). This line is just for testing.
		if (sphere) changeToSphere ();

		double threshold = 0.00001; // 0.01 mm
		
		double aAB,aAP, sAP, sPX, sAX, R, A, lonA2, latA2;
		R = a;
		
		int iterNo = 0;
		int iterates = 1;
		
		
		while (iterates  == 1) {	    	
	    	iterNo = iterNo + 1;
	    	
	    	LineData AB = STX_Inverse (lonA, latA, lonB, latB);
	    	LineData AP = STX_Inverse (lonA, latA, lonP, latP);
	    	
	    	aAB = AB.azi1;
	    	aAP = AP.azi1;
	    	sAP = AP.s12;
	    	A = aAP - aAB;
	    			
	    	sPX = Math.asin(Math.sin(sAP/R)*Math.sin(A))*R;
	    	
	    	if (iterNo == 1)
	    		sAX = 2 * R * Math.atan ( Math.sin((Math.PI/2 + A) / 2) / Math.sin((Math.PI/2-A)/2) * Math.tan((sAP-sPX)/(R*2)));
	    	else
	    		//From second iteration 
	    		sAX = R * Math.atan(Math.cos(A)*Math.tan(sAP/R));
	
	    	LineData AA2 = STX_Direct (lonA, latA, AB.azi1, sAX);
			latA2 = AA2.lat2;
			lonA2 = AA2.lon2;
			
			if (Math.abs(sAX) < threshold) { 
		    				iterates = 0;	
		    				latX = latA2;
		    				lonX = lonA2;
		    		} else {
		    			latA = latA2;
		    			lonA = lonA2;
		    		}   
		    
		}
		
		return iterNo;
	}

}




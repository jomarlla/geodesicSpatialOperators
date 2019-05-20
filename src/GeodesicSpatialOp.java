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
		setAF (6378137, 1 / 298.257223563);
	}
	private double a;
	private double f;
	
	public static double DEG2RAD = Math.PI / 180.0;
	public static double RAD2DEG = 180.0 / Math.PI;
	private boolean localSphere;
	private int nIterations = 1;

	private Geodesic ELLIP;
		    
	private double NA;
	private double NB;
	private double MA;
	private double MB;
	
	public double aAB;
	public double aBA;
	public double sAB;	
	public double latX;
	public double lonX;

	// Parametros en radianes
    public void getRadiusNM (double lonA, double latA, double lonB, double latB) {
    	double exc2 = 2*f-f*f;
    	
    	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
    	aAB = g.azi1 * DEG2RAD;
    	aBA = g.azi2 * DEG2RAD;
    	sAB = g.s12;
    	
    	double slatA = Math.sin(latA); 
    	double slatA2 = slatA * slatA;
   	    double slatB = Math.sin(latB); 
   	    double slatB2 = slatB * slatB;
   	    
   	    NA = a/Math.sqrt(1-exc2*slatA2); // % radius of curvature in the prime vertical
        NB = a/Math.sqrt(1-exc2*slatB2); 
    	
        MA = a*(1-exc2)/ Math.pow ( (1-exc2*slatA2), 1.5); // % radius of the meridian ellipse
        MB = a*(1-exc2)/ Math.pow ( (1-exc2*slatB2), 1.5);
        
        //System.out.println("NA: " + NA + "MA: " + MA);
        //System.out.println("NB: " + NB + "MB: " + MB);
	}
    
    public static double normAzDeg (double az) {
    	az += (az < 0) ? 360 : 0 ;
    	az -= (az > 360) ? 360 : 0 ;
    	
    	return az;
    }
    public static double normAzRad (double az) {
    	double pi2 = 2*Math.PI;
    	az += (az < 0) ? pi2 : 0 ;
    	az -= (az > pi2) ? pi2 : 0 ;
    	
    	return az;
    }
    
    //Unused version, just for checking
    public int[]
    STX_MinDistanceFast2 (double lonA, double latA, double lonB, double latB,
     		            double lonP, double latP) {
    	
   
    	double thresholdlat = 0.000000001 * DEG2RAD; // 1:100000 sec is 0.3 mm
    	double thresholdlon = thresholdlat * Math.cos((latA+latB+latP) / 3.0);
    	
    	int iterNo[] = new int[2];
    	iterNo[0] = 0;
    	int iterates = 1;
    	
    	double lonA0 = lonA;
    	double latA0 = latA;
    	double aAB0 = 0;; 
    	
    	while (iterates  == 1) {
        	double aAB,aAP, sAP, sPX, sAX, R;
        	
        	R = a;
        	iterNo[0]++;
        	
        	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
        	aAB = g.azi1 * DEG2RAD;  
        	if (iterNo[0] == 1) aAB0 = normAzDeg (g.azi1);
  
        	if (aAB < 0) aAB = aAB + 2 * Math.PI;
        	
        	g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latP * RAD2DEG, lonP * RAD2DEG);
        	aAP = g.azi1 * DEG2RAD;
        	sAP = g.s12;    
        	
        	if (aAP < 0) aAP = aAP + 2 * Math.PI;
        	
        	sPX = Math.asin(Math.sin(sAP/R)*Math.sin(aAP-aAB))*R;        	
        	sAX = R*2*Math.atan(Math.sin((Math.PI/2.0+aAP-aAB)/2.0)/Math.sin((Math.PI/2-aAP+aAB)/2.0)*Math.tan((sAP-sPX)/(R*2.0)));
        	//System.out.println("sAx: " + sAX);
        	g = ELLIP.Direct (latA * RAD2DEG, lonA * RAD2DEG, aAB * RAD2DEG, sAX);    		
    		double latA2 = g.lat2 * DEG2RAD;
    		double lonA2 = g.lon2 * DEG2RAD;
    		 
    		
    		if (Math.abs(latA2-latA) < thresholdlat && Math.abs(lonA2-lonA) < thresholdlon) { 
    	    				iterates = 0;	
    	    				latX = latA2;
    	    				lonX = lonA2;
    	    		} else {
    	    			latA = latA2;
    	    			lonA = lonA2;
    	    		}   
    	    
    	}
    	
    	iterates = 1;
    	iterNo[1] = 0;
    	
    	while (iterates  == 1) {
    		iterNo[1]++;
    		
	    	GeodesicData g = ELLIP.Inverse(latX * RAD2DEG, lonX * RAD2DEG, latA0 * RAD2DEG, lonA0 * RAD2DEG);
	    	double dXA = g.s12;
	    	double aXA = normAzDeg (g.azi1);
	    
	    	g = ELLIP.Inverse(latX * RAD2DEG, lonX * RAD2DEG, latP * RAD2DEG, lonP * RAD2DEG);
	    	double sXP = g.s12;
	    	double aXP = normAzDeg (g.azi1);
	    	
	    	double az = 90;
	    	if ( normAzDeg(aXA - aXP) < 180) az = -90;    	
	    	double aXPP = normAzDeg (aXA + az);
	    	
	    	g = ELLIP.Direct (latX * RAD2DEG, lonX * RAD2DEG, aXPP, sXP);
	    	double latX1 = g.lat2;
	    	double lonX1 = g.lon2;
	    	g = ELLIP.Inverse(latP * RAD2DEG, lonP * RAD2DEG, latX1, lonX1);
	    	double dPX1 = g.s12;
	    	
	    	az = normAzDeg (aXA - aXPP);
	    	az += (az > 180) ? -180 : 0;
	    	dPX1 *= (az < 90) ? -1: 1;		
	
	    	g = ELLIP.Direct (latA0 * RAD2DEG, lonA0 * RAD2DEG, aAB0, dXA + dPX1);
	    	lonX = g.lon2 * DEG2RAD;
	    	latX = g.lat2 * DEG2RAD;  
	    	
	    	if (dPX1 < 1E-6 || iterNo[1] > 10) { // less than 1 micron
	    		iterates = 0;
	    	}
    	}
    	
    	return iterNo;
    }
    
    public int
    STX_4PIntersection (double lonA, double latA, double lonB, double latB,
     		            double lonC, double latC,   double lonD, double latD) {
    	
   
    	double thresholdlat = 0.000000001 * DEG2RAD; // 1:100000 sec is 0.3 mm
    	double thresholdlon = thresholdlat * Math.cos((latA+latB+latC+latD) / 4.0);
    	
    	
    	int iterNo = 0;
    	int iterates = 1;
    	
    	double NC, ND, MC, MD, aCD, aDC, sCD;
    	
    	//We never make the calculatios using a local sphere (it is just sperical trigonometry). It is a test.
		if (localSphere) {
			
			double minlon = Math.min (lonA, Math.min(lonB, Math.min(lonC, lonD)));
			double minlat = Math.min (latA, Math.min(latB, Math.min(latC, latD)));
			double maxlon = Math.max (lonA, Math.max(lonB, Math.max(lonC, lonD)));
			double maxlat = Math.max (latA, Math.max(latB, Math.max(latC, lonD)));
			getRadiusNM ( minlon,  minlat,  maxlon,  maxlat);
			this.a = ( NA + NB + MA + MB ) / 4.0;
			this.f = 1 / 1000000000.0;
			
			ELLIP = new Geodesic (a, f);
			
			
			/*
    		getRadiusNM ( lonC,  latC,  lonD,  latD);
    		NC = NA; ND = NB; MC = MA; MD = MB; aCD = aAB; aDC = aBA; sCD = sAB;
    		getRadiusNM ( lonA,  latA,  lonB,  latB);   		
			
			this.a = ( NA + NB + NC + ND + MA + MB + MC + MD ) / 8.0;
			this.f = 1 / 1000000000.0;
			
			ELLIP = new Geodesic (a, f);
			*/
		}
    	
    	while (iterates  == 1) {    		
    		getRadiusNM ( lonC,  latC,  lonD,  latD);
    		NC = NA; ND = NB; MC = MA; MD = MB; aCD = aAB; aDC = aBA; sCD = sAB;
    		getRadiusNM ( lonA,  latA,  lonB,  latB);   		
    		
    		double s2Revision1 = (MC*NA*Math.cos(latA)*NC*Math.cos(latC)*(lonA-lonC)+MA*MC*NC*Math.cos(latC)*Math.tan(aAB)*(latC-latA))/(MC*NA*Math.cos(latA)*Math.sin(aCD)-MA*NC*Math.cos(latC)*Math.cos(aCD)*Math.tan(aAB));
    		double s1, s2;
    		
    		if (aAB == 90 * DEG2RAD || aAB == 270 * DEG2RAD) {
    			s2 = MC/Math.cos(aCD)*(latA-latC);
    			s1 = NA*Math.cos(latA)/Math.sin(aAB)*(lonC-lonA+Math.sin(aCD)/NC/Math.cos(latC)*s2);
    		} else {
    			if (latC == Math.PI * 0.25) {
    				s2 = -(MC*NA*Math.cos(latA)*NC*(lonA-lonC)+MA*MC*NC*Math.tan(aAB)*(latC-latA))/(MA*NC*Math.cos(aCD)*Math.tan(aAB));
    			} else {
    				s2 = s2Revision1;
    			}
    			
    			s1 = MA/Math.cos(aAB)*(latC-latA+Math.cos(aCD)/MC*s2);
    		}
    		
    		iterNo = iterNo + 1;
    		
    		GeodesicData g = ELLIP.Direct(latA * RAD2DEG, lonA * RAD2DEG, aAB * RAD2DEG, s1);
    		double latA2 = g.lat2 * DEG2RAD;
    		double lonA2 = g.lon2 * DEG2RAD;
    		if (latC == 90 * DEG2RAD || latC == -90 * DEG2RAD) {
    			g = ELLIP.Inverse(latC * RAD2DEG, lonC * RAD2DEG, latA2 * RAD2DEG, lonA2 * RAD2DEG);
    			s2 = g.s12;
    		}
    		
    		g = ELLIP.Direct(latC * RAD2DEG, lonC * RAD2DEG, aCD * RAD2DEG, s2);    		
    		double latC2 = g.lat2 * DEG2RAD;
    		double lonC2 = g.lon2 * DEG2RAD;
    		
    		if (Math.abs(latA2-latA) < thresholdlat &&
    		   Math.abs(latC2-latC) < thresholdlat &&
   	           Math.abs(lonA2-lonA) < thresholdlon &&
   	           Math.abs(lonC2-lonC) < thresholdlon) {
    				iterates = 0;	
    				latX = latA2;
    				lonX = lonA2;
    		} else {
    			latA = latA2;
    			lonA = lonA2;
    			latC = latC2;
    			lonC = lonC2;
    		}
    	
    	
    	} //while
    	
    	return iterNo;
    }
    private  void setSpheroid (String args[], int nArgs)  {
    	nIterations = 1;
    	localSphere = false;
    	
    	// Get the spheroid    	
    	if (args[nArgs].equalsIgnoreCase("wgs84")) {
    		a = 6378137;
    		f = 1 / 298.257223563;
    	} else if (args[nArgs].equalsIgnoreCase("grs80")){
    		a = 6378137;
    		f = 1 / 298.257222101;    		
    	}  else if (args[nArgs].equalsIgnoreCase("localsphere")){
    		localSphere = true;
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
    			System.out.println("java -jar GeodesicSpatialOp.jar " + tipo + " loA laA loB laB loP laP esferoid=wgs84,grs80,localsphere [nIterations=1]");
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
    			
    			//Case 3
    			
    			args[0] = tipo;
    			args[1] = "29"; args[2] = "42"; args[3] = "-70"; args[4] = "-35";
    			args[5] = "-22"; args[6] = "64"; 
    			args[7] = "wgs84"; args[8] = "1000";
    			
    			
    			//Case 2
    			/*    			
    			args[0] = tipo;
    			args[1] = "29"; args[2] = "42"; args[3] = "-77"; args[4] = "39";
    			args[5] = "-22"; args[6] = "64"; 
    			args[7] = "wgs84"; args[8] = "1000";
    			*/    			
    			
    			// Case 1
    			/*
    			args[0] = tipo;
    			args[1] = "5"; args[2] = "52"; args[3] = "6"; args[4] = "51.4";
    			args[5] = "5.5"; args[6] = "52"; 
    			args[7] = "wgs84"; args[8] = "1";
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
        	
        	int tipoInt = 3; //MinAcc2 by default
        	
        	if (tipo.equalsIgnoreCase("min")) tipoInt = 1;
        	else if (tipo.equalsIgnoreCase("minAcc")) tipoInt = 2;
        	else tipoInt = 3; //minAcc2
       	
        	int iterNoFast2[] = null;
        	
        	//STX_MinDistanceFast2 and STX_MinDistanceAccurate are not used anymore, because
        	//the new STX_MinDistanceFast gives a high accuracy
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		if (tipoInt == 1) iterNo = jc.STX_MinDistanceFast(xy[0],xy[1],xy[2],xy[3],xy[4],xy[5]);
        		else if (tipoInt == 2) iterNo = jc.STX_MinDistanceAccurate(xy[0],xy[1],xy[2],xy[3],xy[4],xy[5]);
        		else {
        			iterNoFast2 = jc.STX_MinDistanceFast2(xy[0],xy[1],xy[2],xy[3],xy[4],xy[5]);
        			iterNo = iterNoFast2[0];
        		}
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("The algorithm converges in " + iterNo + " iterations");
        	if (tipoInt == 3) {
            	System.out.println("The algorithm converges in " + iterNoFast2[1] + " iterations (improving method)");        		
        	}
        	System.out.println("DDMMSS           LonX: " + ddmmss(jc.lonX * RAD2DEG) + "  LatX: " + ddmmss(jc.latX * RAD2DEG));
        	System.out.println("Radians          LonX: " + (jc.lonX) + "  LatX: " + (jc.latX));
        	System.out.println("Decimal Degrees  LonX: " + (jc.lonX * RAD2DEG) + "  LatX: " + (jc.latX * RAD2DEG));
        	
        	/*
        	double az1,az2;
    	    jc.STX_Inverse(xy[0],xy[1], jc.lonX, jc.latX);
    	    az1 =  jc.aAB;
    	    jc.STX_Inverse(xy[0],xy[1], xy[2],xy[3]);
    	    System.out.println ("x2: " + xy[2] * RAD2DEG);
    	    System.out.println ("x3: " + xy[3] * RAD2DEG);
    	    az2 =  jc.aAB;
    	    System.out.println ("aAB: " + az1 * RAD2DEG);
    	    System.out.println ("aAX: " + az2 * RAD2DEG);
    	    System.out.println ("difaz: " + (az1-az2));
    	    */
    	    
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
        	
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		jc.STX_Direct(xy[0],xy[1],xy[2],xy[3]);
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	
        	double aBA = nAz_rad (jc.aBA);
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("DDMMSS           LonX: " + ddmmss(jc.lonX * RAD2DEG) + "  LatX: " + ddmmss(jc.latX * RAD2DEG) +  " azBA: " + ddmmss(aBA * RAD2DEG));
        	System.out.println("Radians          LonX: " + (jc.lonX) + "  LatX: " + (jc.latX) + "  azBA: " + (aBA));
        	System.out.println("Decimal Degrees  LonX: " + (jc.lonX * RAD2DEG) + "  LatX: " + (jc.latX * RAD2DEG) + "  azBA: " + (aBA * RAD2DEG));
        	
    	
    }
    
    private static double nAz_rad (double az) {
    	double pi2 = 2 * Math.PI;
    	
    	if (az < 0) az += pi2;
    	else if (az > pi2) az -= pi2;
    	
    	return az;
    }
    
    private static double nAz_deg (double az) {
    	double pi2 = 360;
    	
    	if (az < 0) az += pi2;
    	else if (az > pi2) az -= pi2;
    	
    	return az;
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
        	
        	for (int i = 0; i < jc.nIterations; i++) { // 10.000 calculations per second (Core i7 4770)
        		jc.STX_Inverse(xy[0],xy[1],xy[2],xy[3]);
        	}
        	    	    	
        	double runningTime = new Date().getTime() - t0;
        	
        	double aAB = nAz_rad (jc.aAB);
        	double aBA = nAz_rad (jc.aBA);
        	
        	System.out.println("Total Running time in mili seconds: " + runningTime);
        	System.out.println("Number of iterations:" + jc.nIterations);
        	System.out.println("");
        	System.out.println("Distance AB: " + jc.sAB);
        	System.out.println("DDMMSS           azAB: " + ddmmss(aAB * RAD2DEG) +  " azBA: " + ddmmss(aBA * RAD2DEG));
        	System.out.println("Radians          azAB: " + (aAB) + "  azBA: " + (aBA));
        	System.out.println("Decimal Degrees  azAB: " + (aAB * RAD2DEG) + "  azBA: " + (aBA * RAD2DEG));
        	
    	
    }
     
    public static void main(String[] args)  {

    	if (args.length > 0) {
    		if (args[0].equalsIgnoreCase("int")) { mainInt (args); return; }    		  
    		  else if (args[0].equalsIgnoreCase("min")) { mainMin (args); return; } 
    		
    		  //else if (args[0].equalsIgnoreCase("minAcc")) { mainMin (args); return; } 
    		  else if (args[0].equalsIgnoreCase("direct")) { mainDirect (args); return; } 
    		  else if (args[0].equalsIgnoreCase("inverse")) { mainInverse (args); return; } 
    	} 

		System.out.println("Examples:");
		System.out.println("java -jar GeodesicSpatialOp.jar int loA laA loB laB loC laC loD laD wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar min loA laA loB laB loP laP wgs84 [nIterations=1]");
		//System.out.println("java -jar GeodesicSpatialOp.jar minAcc loA laA loB laB loP laP wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar inverse loA laA loB laB wgs84 [nIterations=1]");
		System.out.println("java -jar GeodesicSpatialOp.jar direct loA laA azAB sAB wgs84 [nIterations=1]");
		System.out.println("");
		mainMin(args);
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
			args[1] = "29"; args[2] = "42"; args[3] = "-77"; args[4] = "39";
			args[5] = "0"; args[6] = "6"; args[7] = "-22"; args[8] = "64";
			args[9] = "wgs84"; args[10] = "1000";
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
    	

  
    	
    	System.out.println("Using " + (jc.localSphere ? "a local spehere" : "a spheroid (" + jc.a + "," + jc.f + ")"));
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
    	
    	//Checking the results Line AB
    	//Azimut AB
    	GeodesicData g = jc.ELLIP.Inverse(xy[1] * RAD2DEG, xy[0]  * RAD2DEG, xy[3] * RAD2DEG, xy[2] * RAD2DEG);
    	double dAB = g.s12;
    	double azAB = g.azi1;    	    
    	
    	g = jc.ELLIP.Inverse(xy[1] * RAD2DEG, xy[0]  * RAD2DEG, jc.latX * RAD2DEG, jc.lonX * RAD2DEG);
    	double dAX = g.s12;
    	double azAX = g.azi1;
    	
    	g = jc.ELLIP.Direct(xy[1] * RAD2DEG, xy[0] * RAD2DEG, azAB, dAX);
    	double laPX1 = g.lat2;
    	double loPX1 = g.lon2;
    	
    	System.out.println("");
    	System.out.println("Checking the geodesic AB");
    	System.out.println("Azimuth AB: " + azAB);
    	System.out.println("Azimuth AX: " + azAX);
    	System.out.println("Distance AX: " + dAX);
    	System.out.println("Direct Problem (A, azAB, sAX)  LonPX1: " + loPX1 + "   LatPX1: " + laPX1);
    	
    	//Checking the results Line CD
    	//Azimut CD
    	g = jc.ELLIP.Inverse(xy[5] * RAD2DEG, xy[4]  * RAD2DEG, xy[7] * RAD2DEG, xy[6] * RAD2DEG);
    	double dCD = g.s12;
    	double azCD = g.azi1;    	    
    	
    	g = jc.ELLIP.Inverse(xy[5] * RAD2DEG, xy[4]  * RAD2DEG, jc.latX * RAD2DEG, jc.lonX * RAD2DEG);
    	double dCX = g.s12;
    	double azCX = g.azi1;
    	
    	g = jc.ELLIP.Direct(xy[5] * RAD2DEG, xy[4] * RAD2DEG, azCD, dCX);
    	double laPX2 = g.lat2;
    	double loPX2 = g.lon2;
    	
    	System.out.println("");
    	System.out.println("Checking the geodesic CD");
    	System.out.println("Azimuth CD: " + azCD);
    	System.out.println("Azimuth CX: " + azCX);
    	System.out.println("Distance CX: " + dCX);
    	System.out.println("Direct Problem (C, azCD, sCX)  LonPX2: " + loPX2 + "   LatPX2: " + laPX2);
    	
   	
    	
    	//Distance PX1PX"
    	g = jc.ELLIP.Inverse(laPX1, loPX1, laPX2, loPX2);
    	double dPX1PX2 = g.s12;
    	System.out.println("Distance PX1PX2: " + dPX1PX2);
    	
    }
    
    public static String ddmmss (double g) {
  
    	String ddmmss;
    	int d, m;
    	double t1, s;
    	
    	d = (int) g;  // Truncate the decimals
    	t1 = (g - d) * 60;
    	m = (int) t1;
    	s = (t1 - m) * 60;
    	
    	ddmmss = "" + d + " " + m + " " + String.format("%.8f", s);
    	
    	return ddmmss;
    	
    }
    
    public void
    STX_Inverse (double lonA, double latA, double lonB, double latB) {    	
    	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
    	aAB = g.azi1 * DEG2RAD;
    	aBA = g.azi2 * DEG2RAD;
    	sAB = g.s12;    
    }
    
    public void
    STX_Direct (double lonA, double latA, double azAB, double sAB) {    
    	GeodesicData g = ELLIP.Direct (latA * RAD2DEG, lonA * RAD2DEG, azAB * RAD2DEG, sAB);    		
		latX = g.lat2 * DEG2RAD;
		lonX = g.lon2 * DEG2RAD;   
		aBA = g.azi2 * DEG2RAD; 
    }
    
    //Unused version, just for checking
    public int
    STX_MinDistanceAccurate (double lonA, double latA, double lonB, double latB,
     		            double lonP, double latP) {
    	
   
    	double threshold = 0.000000001 * DEG2RAD; // 1:100000 sec is 0.3 mm
    	double az, aPX, sPX;

    	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
    	aAB = g.azi1 * DEG2RAD; aAB += (aAB < 0) ? Math.PI*2.0 : 0;
    	aBA = (g.azi2 - 180) * DEG2RAD; aBA += (aBA < 0) ? Math.PI*2.0 : 0;
    	sAB = g.s12;    
    	
    	if (aAB  < Math.PI) {
    		   az=(aAB+(aBA-Math.PI))/2;
    		   aPX=az+3*Math.PI/2;
    	} else {
    		   az=(aBA+(aAB-Math.PI))/2;
    		    aPX=az+Math.PI/2;
    	}

    	sPX = 10000;
    	
    	
    	int iterNo = 0;
    	int iterates = 1;
    	
    	double aXP, aXA, aXB, sXP, corrAz;
    	
    	while (iterates  == 1) {
    		iterNo++;
    		
        	g = ELLIP.Direct (latP * RAD2DEG, lonP * RAD2DEG, aPX * RAD2DEG, sPX);    		
    		double latD = g.lat2 * DEG2RAD;
    		double lonD = g.lon2 * DEG2RAD;    		
    		
    		STX_4PIntersection (lonA, latA, lonB, latB, lonP, latP, lonD, latD);
    		
        	g = ELLIP.Inverse(latX * RAD2DEG, lonX * RAD2DEG, latP * RAD2DEG, lonP * RAD2DEG);
        	sXP = g.s12;
        	aXP = g.azi1 * DEG2RAD;
        	aPX = g.azi2 * DEG2RAD; aPX += (aPX < 0) ? Math.PI : -Math.PI;
        	  
        	g = ELLIP.Inverse(latX * RAD2DEG, lonX * RAD2DEG, latA * RAD2DEG, lonA * RAD2DEG);
        	aXA = g.azi1 * DEG2RAD;
        	corrAz=-(aXP-aXA-Math.round((aXP-aXA)/(Math.PI/2))*Math.PI/2);
        	
        	/*
        	  g = ELLIP.Inverse(latX * RAD2DEG, lonX * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
        	  aXB = g.azi1 * DEG2RAD;
        	  corrAz=-(aXP-aXB-Math.round((aXP-aXB)/(Math.PI/2))*Math.PI/2);        	 
        	*/
        	 
            aPX=aPX+corrAz;
        	
            if (Math.abs(corrAz) < threshold || iterNo > 100) iterates = 0;
            
    		sPX = sXP;   	    		       	    
    	}
    	
    	
    	return iterNo;
    }

	// Parametros en radianes
	public int
	STX_MinDistanceFast (double lonA, double latA, double lonB, double latB,
	 		            double lonP, double latP) {
		
	
		double thresholdlat = 0.000000001 * DEG2RAD; // 1:100000 sec is 0.3 mm
		double thresholdlon = thresholdlat * Math.cos((latA+latB+latP) / 3.0);
		
		
		int iterNo = 0;
		int iterates = 1;
		
		while (iterates  == 1) {
	    	double aAB,aAP, sAP, sPX, sAX, R;
	    	
	    	R = a;
	    	iterNo = iterNo + 1;
	    	
	    	GeodesicData g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latB * RAD2DEG, lonB * RAD2DEG);
	    	aAB = g.azi1 * DEG2RAD;  	
	
	    	if (aAB < 0) aAB = aAB + 2 * Math.PI;
	    	
	    	g = ELLIP.Inverse(latA * RAD2DEG, lonA * RAD2DEG, latP * RAD2DEG, lonP * RAD2DEG);
	    	aAP = g.azi1 * DEG2RAD;
	    	sAP = g.s12;    
	    	
	    	if (aAP < 0) aAP = aAP + 2 * Math.PI;
	    	
	    	sPX = Math.asin(Math.sin(sAP/R)*Math.sin(aAP-aAB))*R;        	
	    	
	    	if (iterNo == 1)
	    		sAX = R*2*Math.atan(Math.sin((Math.PI/2.0+aAP-aAB)/2.0)/Math.sin((Math.PI/2-aAP+aAB)/2.0)*Math.tan((sAP-sPX)/(R*2.0)));
	    	else
	    		//From second iteration 
	    		sAX = R*Math.atan(Math.cos(aAP-aAB)*Math.tan(sAP/R));

	    	g = ELLIP.Direct (latA * RAD2DEG, lonA * RAD2DEG, aAB * RAD2DEG, sAX);    		
			double latA2 = g.lat2 * DEG2RAD;
			double lonA2 = g.lon2 * DEG2RAD;
			
			
			if (Math.abs(latA2-latA) < thresholdlat && Math.abs(lonA2-lonA) < thresholdlon) { 
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




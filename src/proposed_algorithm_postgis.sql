/*
 * PlpgSQL for PostgreSQL implementation of the intersection and point-to-line solutions 
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

/* 
  Get the point X defined by the intersection of the geodesic pa-pb with the geodesic pc-pd
  geography STX_GeodesicIntersection ( pa geography, pb geography, pc geography, pd geography) RETURNS Geography


  Get the point on the geodesic line pa-pb closest to the point pp
  geography STX_GeodesicMinDistance (pa geography, pb geography, pp geography) RETURNS Geography AS
  
  The distance itself will be: ST_Distance (pp, STX_GeodesicMinDistance (pa geography, pb geography, pp geography) )
  
  Private functions (do not use):
  _getRadiusNM     
*/

BEGIN;

CREATE OR REPLACE FUNCTION _getRadiusNM (a float8, f float8, lonA float8, latA float8, lonB float8, latB float8,
  out NA float8, out NB float8, out MA float8, out MB float8, out aAB float8, out aBA float8, out sAB float8, out pa geography, out pb geography
) RETURNS record AS
$$
DECLARE
  exc2 float8;
  
  slatA float8;
  slatA2 float8;
  slatB float8;
  slatB2 float8;
  
BEGIN

	  exc2:= 2*f-f*f;

	  pa := ST_MakePoint (degrees(lonA), degrees(latA))::geography;
	  pb := ST_MakePoint (degrees(lonB), degrees(latB))::geography;
	  
	  --RAISE NOTICE 'X (%)', ST_X(pa::geometry);
	  
	  sAB:= ST_Distance (pa,pb); 
	  aAB:= ST_Azimuth (pa,pb);
	  aBA:= ST_Azimuth (pb,pa); --returns radians
	  
	  IF aBA > radians(180) THEN
	    aBA:= aBA - radians(180);
	  ELSE  
	    aBA:= aBA + radians(180);		
	  END IF;
    
	  slatA := sin(latA); slatA2 = slatA * slatA;
	  slatB := sin(latB); slatB2 = slatB * slatB;
	  
	  NA:= a/sqrt(1-exc2*slatA2); -- % radius of curvature in the prime vertical
      NB:= a/sqrt(1-exc2*slatB2);      
	  MA:= a*(1-exc2)/ power ( (1-exc2*slatA2), 1.5); -- % radius of the meridian ellipse
      MB:= a*(1-exc2)/ power ( (1-exc2*slatB2), 1.5);
END;
$$ LANGUAGE plpgsql IMMUTABLE STRICT;
      


CREATE OR REPLACE FUNCTION STX_GeodesicIntersection (
  pa geography, pb geography, pc geography, pd geography
) RETURNS Geography AS
$$
DECLARE

  a float8;
  f float8;
  thresholdlat float8;
  thresholdlon float8;
  latA float8;
  lonA float8;
  latB float8;
  lonB float8;
  latC float8;
  lonC float8;
  latD float8;
  lonD float8;
  
  pa2 geometry;
  pc2 geometry;
  latA2 float8;
  lonA2 float8;
  latC2 float8;
  lonC2 float8;
  
  sAB float8;
  aAB float8;
  aBA float8;
  NA float8;
  NB float8;
  MA float8;
  MB float8;

  sCD float8;
  aCD float8;
  aDC float8;
  NC float8;
  ND float8;
  MC float8;
  MD float8;
  
  iterNo integer;
  iterates integer;
  
  s1 float8;
  s2 float8;
  
  rec1 record;
  rec2 record;
  s2Revision1 float8;
  
BEGIN
	lonA:= radians(st_x(pa::geometry)); latA:= radians(st_y(pa::geometry));
	lonB:= radians(st_x(pb::geometry)); latB:= radians(st_y(pb::geometry));
	lonC:= radians(st_x(pc::geometry)); latC:= radians(st_y(pc::geometry));
	lonD:= radians(st_x(pd::geometry)); latD:= radians(st_y(pd::geometry));
	
	a:= 6378137;
	f:= 1/298.257223563;  -- % for GRS80 use f=1/298.257222100882711243
	thresholdlat:= 0.00001/206265; -- % 0.00001" corresponds to 0.3 mm
	thresholdlon:= thresholdlat*cos((latA+latB+latC+latD) / 4.0);
	
	iterNo:= 0;
    iterates:= 1;
	
	WHILE iterates = 1 LOOP	
	  rec1:= _getRadiusNM (a, f, lonA, latA, lonB, latB);
	  NA:= rec1.NA; NB:= rec1.NB;
	  MA:= rec1.MA; MB:= rec1.MB;
	  sAB:= rec1.sAB; aAB:= rec1.aAB; aBA:= rec1.aBA;
	  pa:= rec1.pa; pb:= rec1.pb;
	  
	  rec2:= _getRadiusNM (a, f, lonC, latC, lonD, latD);
	  NC:= rec2.NA; ND:= rec2.NB;
	  MC:= rec2.MA; MD:= rec2.MB;
	  sCD:= rec2.sAB; aCD:= rec2.aAB; aDC:= rec2.aBA;
	  pc:= rec2.pa; pd:= rec2.pb;
	  
	  --RAISE NOTICE 'RADIUS: % % % % % % %',NA,NB,MA,MB,sAB, aAB, aBA;
	  --RAISE NOTICE 'RADIUS: % % % % % % %',NC,ND,MC,Md,sCD, aCD, aDC;
	  
	  s2Revision1:= (MC*NA*cos(latA)*NC*cos(latC)*(lonA-lonC)+MA*MC*NC*cos(latC)*tan(aAB)*(latC-latA))/(MC*NA*cos(latA)*sin(aCD)-MA*NC*cos(latC)*cos(aCD)*tan(aAB));
	  
	  IF aAB = radians(90) OR aAB = radians(270) THEN
	    s2:= MC/cos(aCD)*(latA-latC);
		s1:= NA*cos(latA)/sin(aAB)*(lonC-lonA+sin(aCD)/NC/cos(latC)*s2);
	  ELSE  
	    IF latC = radians(90) THEN
		  s2:= -(MC*NA*cos(latA)*NC*(lonA-lonC)+MA*MC*NC*tan(aAB)*(latC-latA))/(MA*NC*cos(aCD)*tan(aAB));
		ELSE 
		  s2:= s2Revision1;
		END IF;
		
		s1:= MA/cos(aAB)*(latC-latA+cos(aCD)/MC*s2);
	  END IF;
	
	  
      iterNo:= iterNo+1;
	  
	  pa2:= ST_Project (pa, s1, aAB)::geometry;
	  latA2:= radians (ST_Y(pa2));
	  lonA2:= radians (ST_X(pa2));
	  
	  IF latC = radians (90) OR latC = radians (-90) THEN	  
	    s2:= ST_Distance (pc, ST_MakePoint (degrees(lonA2), degrees(latA2))::geography);	  
	  END IF;
	
	  pc2:= ST_Project (pc, s2, aCD)::geometry; 
	  latC2:= radians (ST_Y(pc2));
	  lonC2:= radians (ST_X(pc2));
	
	  
      IF abs(latA2-latA)<thresholdlat AND abs(latC2-latC)<thresholdlat 
	     AND abs(lonA2-lonA)<thresholdlon AND abs(lonC2-lonC)<thresholdlon THEN		 		 
        iterates=0;  -- % end of iterations
      ELSE
        latA=latA2;
        lonA=lonA2;
        latC=latC2;
        lonC=lonC2;
      END IF;
    
 	  
	END LOOP;
	
	--RAISE NOTICE 'ITERATIONS: %', iterNo;
	
	RETURN ST_MakePoint (degrees(lonA2), degrees(latA2))::geography;
	

END;
$$ LANGUAGE plpgsql IMMUTABLE STRICT;


/*
s1=# select st_aslatlontext(stx_intersection ('POINT (29 42)'::geography,'POINT (-77 39)'::geography,'POINT (0 6)'::geography,'POINT (-22 64)'::geography)::geometry,intersection);
       st_aslatlontext
-----------------------------
 54 43 1.3066 -14 33 49.8807
 */
 
/*
s1=# select st_aslatlontext(stx_intersection ('POINT (-92 35)'::geography,'POINT (52 40)'::geography,'POINT (20 -8)'::ge
ography,'POINT (-95 49)'::geography)::geometry,'D M S.SSSS');
       st_aslatlontext
------------------------------
 50 28 44.7508 -79 16 58.0861
*/

CREATE OR REPLACE FUNCTION STX_GeodesicIntersection (
  pa geometry, pb geometry, pc geometry, pd geometry
) RETURNS geometry AS 'select STX_GeodesicIntersection ($1::geography, $2::geography, $3::geography, $4::geography)::geometry'
  LANGUAGE SQL IMMUTABLE STRICT;



CREATE OR REPLACE FUNCTION STX_GeodesicMinDistance (
  pa geography, pb geography, pp geography
) RETURNS Geography AS
$$
DECLARE
    a float8;
    f float8;
    exc2 float8;

    thresholdlat float8;
    thresholdlon float8;
    latA float8;
    lonA float8;
    latB float8;
    lonB float8;
    latP float8;
    lonP float8;
	latA2 float8;
	lonA2 float8;
	
	R float8;
	
	iterates integer;
	iterNo integer;
	
	sAB float8;
	aAB float8;
	sAP float8;
	aAP float8;
	sAX float8;
	sPX float8;
	
	pa2 geometry;
	
BEGIN
	lonA:= radians(st_x(pa::geometry)); latA:= radians(st_y(pa::geometry));
	lonB:= radians(st_x(pb::geometry)); latB:= radians(st_y(pb::geometry));
	lonP:= radians(st_x(pp::geometry)); latP:= radians(st_y(pp::geometry));

	a:= 6378137;
	f:= 1/298.257223563;  -- % for GRS80 use f=1/298.257222100882711243
	thresholdlat:= 0.00001/206265; -- % 0.00001" corresponds to 0.3 mm
	thresholdlon:= thresholdlat*cos((latA+latB+latP) / 3.0);
	R:= a;
	
    iterates:= 1;
	iterNo:= 0;
	
  	WHILE iterates = 1 LOOP		
		pa := ST_MakePoint (degrees(lonA), degrees(latA))::geography;

	    -- sAB:= ST_Distance (pa,pb);  --no se utiliza
	    aAB:= ST_Azimuth (pa,pb); --aBA no se utiliza
		IF aAB < 0 THEN
		  aAB:= aAB + 2 * pi();		 
		END IF;
		
        sAP:= ST_Distance (pa,pp); 
	    aAP:= ST_Azimuth (pa,pp); --aPA no se utiliza
		IF aAP < 0 THEN
		  aAP:= aAP + 2 * pi();		 
		END IF;		
		
		sPX:= asin(sin(sAP/R)*sin(aAP-aAB))*R;
		
	    IF iterNo = 1 THEN
		  sAX:= R*2*atan(sin((pi()/2.0+aAP-aAB)/2.0)/sin((pi()/2-aAP+aAB)/2.0)*tan((sAP-sPX)/R/2.0));
		ELSE
	      sAX = R*atan(cos(aAP-aAB)*tan(sAP/R));
		END IF;
		
	    pa2:= ST_Project (pa, sAX, aAB)::geometry; 
	    latA2:= radians (ST_Y(pa2));
	    lonA2:= radians (ST_X(pa2));
		
		IF ABS (latA2-latA) < thresholdlat AND ABS (lonA2-lonA) < thresholdlon THEN
          iterNo:= iterNo+1;
          iterates:= 0; -- % end of iterations		  		  
		ELSE
          iterNo:= iterNo+1;
          latA:= latA2;
          lonA:= lonA2;
		END IF;
	END LOOP;
	
	--RAISE NOTICE 'latX: %', degrees(latA);
	--RAISE NOTICE 'ITERATIONS: %', iterNo;

    RETURN ST_MakePoint (degrees(lonA2), degrees(latA2))::geography;


END;
$$ LANGUAGE plpgsql IMMUTABLE STRICT;
	 
/*
s1=# select st_aslatlontext(stx_3pdistance ('POINT (5 52)'::geography,'POINT (6 51.4)'::geography,'POINT (5.5 52)'::geography)::geometry,'D M S.SSSS');
      st_aslatlontext
----------------------------
 51 50 45.9212 5 15 37.5426
 
s1=# select st_aslatlontext(stx_3pdistance ('POINT (29 42)'::geography,'POINT (-77 39)'::geography,'POINT (-22 64)'::geography)::geometry,'D M S.SSSSS');
        st_aslatlontext
--------------------------------
 54 55 42.71339 -21 56 14.24782
 
		 s1=# select (st_azimuth ('POINT (29 42)'::geography,'POINT (-77 39)'::geography));
			 st_azimuth
		--------------------
		 -0.884772900760781
		 
		 s1=# select st_azimuth ('POINT(29 42)'::geography, px) from (select stx_3pdistance ('POINT (29 42)'::geography,'POINT (-
		77 39)'::geography,'POINT (-22 64)'::geography) as px) as tabla;
			st_azimuth
		-------------------
		 -0.88477290076078
		 
			s1=# select st_aslatlontext (st_makepoint (180+degrees(-0.88477290076078),0)::geometry,'D M S.SSSSS');
				   st_aslatlontext
			-----------------------------
			 0 0 0.00000 129 18 22.48905

*/ 
	  

	    
CREATE OR REPLACE FUNCTION STX_GeodesicMinDistance (
  pa geometry, pb geometry, pp geometry
) RETURNS geometry AS 'select STX_GeodesicMinDistance ($1::geography, $2::geography, $3::geography)::geometry'
  LANGUAGE SQL IMMUTABLE STRICT;

END;

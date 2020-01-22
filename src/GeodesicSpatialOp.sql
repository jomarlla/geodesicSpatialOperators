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
*/

BEGIN;


CREATE OR REPLACE FUNCTION STX_GeodesicIntersection (
  pa geography, pb geography, pc geography, pd geography
) RETURNS Geography AS
$$
DECLARE
  R float8;
  threshold float8;
  latA float8;
  lonA float8;
  latC float8;
  lonC float8;
  latB float8;
  latD float8;
  
  pa2 geometry;
  pc2 geometry;
  latA2 float8;
  lonA2 float8;
  latC2 float8;
  lonC2 float8;
  
  sAB float8;
  sCD float8;
  sAC float8;
  sAX float8;
  sCX float8;

  aAB float8;
  aCD float8;
  aAC float8;
  aCA float8;
  
  
  A float8;
  C float8;
  angle1 float8;
  angle2 float8;
  
  iterNo integer;
  iterates integer;
  
BEGIN
	lonA:= radians(st_x(pa::geometry)); latA:= radians(st_y(pa::geometry));
	lonC:= radians(st_x(pc::geometry)); latC:= radians(st_y(pc::geometry));
	latB:= radians(st_y(pb::geometry)); latD:= radians(st_y(pd::geometry));
	
	R:= 6378137;
	threshold:= 0.0001; -- 0.1 mm
	
	iterNo:= 0;
    iterates:= 1;
	
	WHILE iterates = 1 LOOP	
	
	  aAB:= ST_Azimuth (pa,pb);	  
	  aCD:= ST_Azimuth (pc,pd);
	  aAC:= ST_Azimuth (pa,pc); -- These three should be just one sentence
	  sAC:= ST_Distance (pa,pc); 
	  aCA:= ST_Azimuth (pc,pa); 
	  
	  A:= aAC - aAB;	  	  
	  C:= aCD - aCA; 
	  
      angle1:= atan(sin(sAC/R)/(cos(aAC/R)*cos(A)+(1/tan(C))*sin(A)));
	  sAX:= R * angle1;
	  
	  angle2:= atan(sin(sAC/R)/(cos(aAC/R)*cos(C)+(1/tan(A))*sin(C)));
	  sCX:= R * angle2;
	  
	  iterNo:= iterNo+1;
	  
	  pa2:= ST_Project (pa, sAX, aAB)::geometry;
	  latA2:= radians (ST_Y(pa2));
	  lonA2:= radians (ST_X(pa2));
	  
	  -- RAISE NOTICE 'Iteration: %', iterNo;
	  -- RAISE NOTICE 'sAX: %', sAX;
	  -- RAISE NOTICE 'sCX: %', sCX;
 	  -- RAISE NOTICE 'latA2: %', ST_Y(pa2);
	  -- RAISE NOTICE 'lonA2: %', ST_X(pa2);
	  
	  pc2:= ST_Project (pc, sCX, aCD)::geometry;
	  latC2:= radians (ST_Y(pc2));
	  lonC2:= radians (ST_X(pc2));
 	  --RAISE NOTICE 'latC2: %', ST_Y(pc2);
	  --RAISE NOTICE 'lonC2: %', ST_X(pc2);
	  
      pa:= ST_MakePoint (degrees(lonA2), degrees(latA2))::geography;
	  pc:= ST_MakePoint (degrees(lonC2), degrees(latC2))::geography;
		
      IF (abs(sAX)<threshold AND abs(sCX)<threshold) OR
		 (abs(latA2-latC2) < 1E-14 AND abs(lonA2-lonC2) < 1E-14) THEN		 		 
        iterates = 0;  -- % end of iterations
      ELSE
        latA=latA2;
        lonA=lonA2;
        latC=latC2;
        lonC=lonC2;		
      END IF;
    
 	  
	END LOOP;
	
	RAISE NOTICE 'LONX %   LATX %   NITERATIONS: %', degrees(lonA), degrees(latA), iterNo;
	
	RETURN pa;
	

END;
$$ LANGUAGE plpgsql IMMUTABLE STRICT;


/*
WITH Geodesics (ncase, pa, pb, pc, pd) AS 	
   ( SELECT 'Case 1', 'POINT (14.5 54)','POINT (14.6 54.2)', 'POINT (14.4 54.1)','POINT (14.7 54)' union
	 SELECT 'Case 2', 'POINT (5 52)','POINT (6 51.4)', 'POINT (4.5 51.5)','POINT (5.5 52)' union
	 SELECT 'Case 3', 'POINT (29 42)','POINT (-77 39)', 'POINT (0 6)','POINT (-22 64)' union
     SELECT 'Case 4', 'POINT (-92 35)','POINT (52 40)', 'POINT (20 -8)','POINT (-95 49)' union
	 
     SELECT 'Case 5', 'POINT (5 90)','POINT (5 0)', 'POINT (-30 70)','POINT (40 70)' union
     SELECT 'Case 6', 'POINT (-175 80)','POINT (5 0)', 'POINT (-30 60)','POINT (40 80)' union
     SELECT 'Case 7', 'POINT (-170 85)','POINT (12 -15)', 'POINT (-58 26)','POINT (120 75)' union
     SELECT 'Case 8', 'POINT (105 63)','POINT (79 42)', 'POINT (-167 38)','POINT (-100 23)' union
     SELECT 'Case 9', 'POINT (-42 40)','POINT (63 65.5)', 'POINT (-41.8 40)','POINT (62.9 65.6)' ),
    
	Intersections (ncase, px, pa, pb, pc, pd) AS
   ( SELECT ncase, STX_GeodesicIntersection (pa, pb, pc, pd), 
	               pa::geography, pb::geography, pc::geography, pd::geography
     FROM Geodesics )
	 
SELECT ncase, ST_Azimuth (pa, pb) - ST_Azimuth (pa, px) as "aAB-aAX", 
              ST_Azimuth (pc, pd) - ST_Azimuth (pc, px) as "aCD-aCX", 
			  ST_X(px::geometry) as "lonX",
			  ST_Y(px::geometry) as "latY",
			  ST_AsLatLonText ( px, 'D°M''S.SSSSSSSSS"' ) as "DD MM SS"
FROM Intersections ORDER BY ncase;


 ncase  |        aAB-aAX        |        aCD-aCX        |       lonX        |       latY       |                 DD MM SS
--------+-----------------------+-----------------------+-------------------+------------------+-------------------------------------------
 Case 1 |   1.0680345496894e-13 |  4.48530101948563e-14 |  14.5285478498717 | 54.0573013091912 | 54°3'26.284713088" 14°31'42.772259538"
 Case 2 |  1.68753899743024e-14 |  6.66133814775094e-16 |  5.22745711452158 | 51.8656654013763 | 51°51'56.395444955" 5°13'38.845612278"
 Case 3 | -3.33066907387547e-16 |                     0 | -14.5638557443078 | 54.7170296089477 | 54°43'1.306592212" -14°33'49.880679508"
 Case 4 | -1.11022302462516e-16 |                     0 |  -79.282801686624 | 50.4790974467667 | 50°28'44.750808360" -79°16'58.086071846"
 Case 5 |                     0 |  1.77635683940025e-15 |                 5 | 73.4002981961227 | 73°24'1.073506042" 5°0'0.000000000"
 Case 6 |  4.39478391551296e-18 |  -4.9960036108132e-16 |  4.99999999999998 | 77.5238059595557 | 77°31'25.701454400" 4°59'60.000000000"
 Case 7 |   -1.346145417358e-15 |  2.08166817117217e-16 |  34.2479946167131 | 89.5261957615611 | 89°31'34.304741620" 34°14'52.780620167"
 Case 8 |     -3.14159265358979 | -1.33226762955019e-15 | -135.325739600252 | 36.2632698564387 | 36°15'47.771483179" -135°19'32.662560908"
 Case 9 |  1.66533453693773e-16 | -2.22044604925031e-16 |  15.4380985995178 | 68.6469095374995 | 68°38'48.874334998" 15°26'17.154958264"
 
 
 -- OR USING A TABLE AS:
CREATE TABLE Geodesics (ncase integer, pa geography, pb geography, pc geography, pd geography);
INSERT INTO Geodesics VALUES (1, 'POINT(14.5 54)','POINT(14.6 54.2)', 'POINT(14.4 54.1)','POINT(14.7 54)');
INSERT INTO Geodesics VALUES (2, 'POINT(5 52)','POINT(6 51.4)', 'POINT(4.5 51.5)','POINT(5.5 52)');
INSERT INTO Geodesics VALUES (3, 'POINT(29 42)','POINT(-77 39)', 'POINT(0 6)','POINT(-22 64)');
INSERT INTO Geodesics VALUES (4, 'POINT(-92 35)','POINT(52 40)', 'POINT(20 -8)','POINT(-95 49)');

WITH Intersections (ncase, px, pa, pb, pc, pd) AS
   ( SELECT ncase, STX_GeodesicIntersection (pa, pb, pc, pd), pa, pb, pc, pd
     FROM Geodesics )
	 
SELECT ncase, ST_Azimuth (pa, pb) - ST_Azimuth (pa, px) as "aAB-aAX", 
              ST_Azimuth (pc, pd) - ST_Azimuth (pc, px) as "aCD-aCX", 
			  ST_X(px::geometry) as "lonX",
			  ST_Y(px::geometry) as "latY",
			  ST_AsLatLonText ( px::geometry, 'D°M''S.SSSSSSSSS"' ) as "DD MM SS"
FROM Intersections ORDER BY ncase;
*/


/*
SELECT ST_AsLatLonText ( STX_GeodesicIntersection 
	('POINT (29 42)','POINT (-77 39)', 'POINT (0 6)','POINT (-22 64)'), 'D°M''S.SSSSSSSSS"' );

	             st_aslatlontext
-----------------------------------------
 54°43'1.306592212" -14°33'49.880679508"
*/
 


CREATE OR REPLACE FUNCTION STX_GeodesicIntersection (
  pa geometry, pb geometry, pc geometry, pd geometry
) RETURNS geometry AS 'select STX_GeodesicIntersection ($1::geography, $2::geography, $3::geography, $4::geography)::geometry'
  LANGUAGE SQL IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION STX_GeodesicIntersection (
  pa text, pb text, pc text, pd text
) RETURNS geometry AS 'select STX_GeodesicIntersection ($1::geography, $2::geography, $3::geography, $4::geography)::geometry'
  LANGUAGE SQL IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION STX_GeodesicMinDistance (
  pa geography, pb geography, pp geography
) RETURNS Geography AS
$$
DECLARE
    f float8;
    exc2 float8;

    threshold float8;
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
	A float8;
	
	pa2 geometry;
	
BEGIN
	lonA:= radians(st_x(pa::geometry)); latA:= radians(st_y(pa::geometry));
	lonB:= radians(st_x(pb::geometry)); latB:= radians(st_y(pb::geometry));
	lonP:= radians(st_x(pp::geometry)); latP:= radians(st_y(pp::geometry));

	R:= 6378137;
	threshold:= 0.0001; -- 0.1 mm	
	
    iterates:= 1;
	iterNo:= 0;
	
  	WHILE iterates = 1 LOOP		
		iterNo:= iterNo+1;
		
		pa := ST_MakePoint (degrees(lonA), degrees(latA))::geography;

	    aAB:= ST_Azimuth (pa,pb); --aBA not needed
        sAP:= ST_Distance (pa,pp); 
	    aAP:= ST_Azimuth (pa,pp); --aPA not needed
		
		A:= aAP - aAB;
		
		sPX:= asin(sin(sAP/R)*sin(A))*R;
						
	    IF iterNo = 1 THEN
		  sAX:= R*2*atan(sin((pi()/2 + A) / 2) / sin((pi()/2 - A) / 2) * tan((sAP-sPX)/(2*R)));
		ELSE
	      sAX = R*atan(cos(A)*tan(sAP/R));
		END IF;
		
	    pa2:= ST_Project (pa, sAX, aAB)::geometry; 
	    latA2:= radians (ST_Y(pa2));
	    lonA2:= radians (ST_X(pa2));
		
		IF ABS (sAX) < threshold THEN
          iterates:= 0; -- % end of iterations		  		  
		ELSE
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
WITH PX (geom) AS (
  SELECT STX_GeodesicMinDistance ('POINT (29 42)','POINT (-70 -35)', 'POINT (-22 64)')
) SELECT ST_AsLatLonText (geom, 'D°M''S.SSSSSSSSS"') as PX,
    ST_Distance (geom, 'POINT (-22 64)'::geography) as dXP,
    degrees (ST_Azimuth (geom, 'POINT (29 42)'::geography) -
		     ST_Azimuth (geom, 'POINT (-22 64)'::geography) ) as "aXA-aXP"
  FROM PX;
*/

/*
s1=# select st_aslatlontext(STX_GeodesicMinDistance ('POINT (5 52)'::geography,'POINT (6 51.4)'::geography,'POINT (5.5 52)'::geography)::geometry,'D M S.SSSS');
      st_aslatlontext
----------------------------
 51 50 45.9212 5 15 37.5426
 
s1=# select st_aslatlontext(STX_GeodesicMinDistance ('POINT (29 42)'::geography,'POINT (-77 39)'::geography,'POINT (-22 64)'::geography)::geometry,'D M S.SSSSS');
        st_aslatlontext
--------------------------------
 54 55 42.71339 -21 56 14.24782
 
		 s1=# select (st_azimuth ('POINT (29 42)'::geography,'POINT (-77 39)'::geography));
			 st_azimuth
		--------------------
		 -0.884772900760781
		 
		 s1=# select st_azimuth ('POINT(29 42)'::geography, px) from (select STX_GeodesicMinDistance ('POINT (29 42)'::geography,'POINT (-
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
  
CREATE OR REPLACE FUNCTION STX_GeodesicMinDistance (
  pa text, pb text, pp text
) RETURNS geometry AS 'select STX_GeodesicMinDistance ($1::geography, $2::geography, $3::geography)::geometry'
  LANGUAGE SQL IMMUTABLE STRICT;  

END;

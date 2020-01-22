/* Copyright (c) Jose Martinez-Llario (2018-2019) and licensed
 * under the GPL-3.0-or-later License. 
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 */
	

  /*
     Run first the file "GeodesicSpatialOp.sql" which is the implementation of the presented algorithm designed by 
	 Jose Martinez-Llario and Sergio Baselga.
	 It will add the stored procedure "STX_GeodesicIntersection" that its needed in this code to get the true geodesic intersection.
	 
	 Run the file "table_geodesicintersection_several_software.sql" which created the tables for the data
  */
  
 
-- PX1,PX2,PX3,PX4 with the paper proposed algorithm in PostGIS 
-- The function STX_GeodesicIntersection is a SQL stored procedure implementing the algorithm from the authors. You must run the appropriate sql file to create this function.

insert into inter (ncase, software, geom)
select l1.ncase, 'PG_Paper', 
 STX_GeodesicIntersection (st_startpoint (l1.geom), st_endpoint (l1.geom),
                           st_startpoint (l2.geom), st_endpoint (l2.geom))
from lines l1, lines l2 where  l1.name = 'AB' and l2.name = 'CD' and l1.ncase = l2.ncase;

-- Print out all the X coordinates in DMS and decimal degrees  
		/* -- S.12345 (0.3mm in equador) S.123456789 (0.03um in equador)
		   
		SELECT i.ncase, i.software, st_astext(geom) as degrees, (ST_AsLatLonText(geom, 'D°M''S.SSSSSSSSS"')) as DMS from inter i where i.software = 'PG_Paper' order by i.software, i.ncase;	
		
		 ncase  | software |                  degrees                  |                   dms
		--------+----------+-------------------------------------------+------------------------------------------
		 Case 1 | PG_Paper | POINT(14.5285478498717 54.0573013091912)  | 54°3'26.284713088" 14°31'42.772259538"
		 Case 2 | PG_Paper | POINT(5.22745711452158 51.8656654013763)  | 51°51'56.395444955" 5°13'38.845612278"
		 Case 3 | PG_Paper | POINT(-14.5638557443078 54.7170296089477) | 54°43'1.306592212" -14°33'49.880679508"
		 Case 4 | PG_Paper | POINT(-79.282801686624 50.4790974467667)  | 50°28'44.750808360" -79°16'58.086071846"
		 Case 5 | PG_Paper | POINT(4.99999999999998 73.4002981961227)  | 73°24'1.073506042" 4°59'60.000000000"
				

		*/
/* Example for the paper 

		WITH intersectionTable(ncase, geom) AS (
		  SELECT l1.ncase,
			  STX_GeodesicIntersection (st_startpoint (l1.geom), st_endpoint (l1.geom),
										st_startpoint (l2.geom), st_endpoint (l2.geom))
		  FROM lines l1, lines l2 
		  WHERE l1.name = 'AB' and l2.name = 'CD' and l1.ncase = l2.ncase)
		  
		SELECT ncase, ST_astext(geom) as degrees, 
			   ST_AsLatLonText(geom, 'D°M''S.SSSSSSSSS"') as DMS FROM intersectionTable;
		 ncase  |                  degrees                  |                   dms
		--------+-------------------------------------------+------------------------------------------
		 Case 1 | POINT(14.5285478498717 54.0573013091912)  | 54°3'26.284713088" 14°31'42.772259538"
		 Case 2 | POINT(5.22745711452158 51.8656654013763)  | 51°51'56.395444955" 5°13'38.845612278"
		 Case 3 | POINT(-14.5638557443078 54.7170296089477) | 54°43'1.306592212" -14°33'49.880679508"
		 Case 4 | POINT(-79.282801686624 50.4790974467667)  | 50°28'44.750808360" -79°16'58.086071846"
		 Case 5 | POINT(4.99999999999998 73.4002981961227)  | 73°24'1.073506042" 4°59'60.000000000"
			   
*/
	   
 
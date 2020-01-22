/*
 * Experiments for the paper [1]
 * [1] New paper (sending status)
 *
 * Copyright (c) Jose Martinez-Llario (2018-2019) and licensed
 * under the GPL-3.0-or-later License. 
 * The algorithm uses the library geographiclib from Charles Karney to
 * perform the direct and inverse geodetic problems
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
  Instructions:
  
  1 Create a new PostGIS database
  2 Run the file "GeodesicSpatialOp.sql" in the database to add the new spatial algorithm 
  3 Run this file in the database.

  The table "points" contains the points a,b,p for the 3 test cases
  The table "pointongeodesic" contains the minimum geodesic distance calculated by Oracle, Google Earth, ArcGIS and PostGIS:
  
  	If you want to get more information read these:
		Oracle Spatial: instructions in "geodesicmindistance_with_oracle.sql"
		Google Earth Engine: instructions in "geodesicmindistance_with_google_earthengine.js"
		PostGIS with local projections: from lines 48-50 in this file
		PostGIS with the new proposed algorithm: from lines 52-60
*/

create table points (gid serial primary key, ncase varchar, geoma geometry (point, 4326), geomb geometry (point, 4326), geomp geometry (point, 4326));
insert into points (ncase, geoma, geomb, geomp) values ('Case 1', 'SRID=4326;POINT (5 52)', 'SRID=4326;POINT (6 51.4)', 'SRID=4326;POINT (5.5 52)');
insert into points (ncase, geoma, geomb, geomp) values ('Case 2', 'SRID=4326;POINT (29 42)', 'SRID=4326;POINT (-77 39)', 'SRID=4326;POINT (-22 64)');
insert into points (ncase, geoma, geomb, geomp) values ('Case 3', 'SRID=4326;POINT (29 42)', 'SRID=4326;POINT (-70 -35)', 'SRID=4326;POINT (-22 64)');
--Extra cases
insert into points (ncase, geoma, geomb, geomp) values ('Case 4', 'SRID=4326;POINT (29 42)', 'SRID=4326;POINT (-70 -35)', 'SRID=4326;POINT (-22 64)');

create table pointongeodesic (gid serial primary key, ncase varchar, software varchar, distance double precision, geom geometry (point, 4326));

-- PostGIS transformation to local projection
insert into pointongeodesic (ncase, software, distance, geom) 
  select ncase,'PostGIS', ST_Distance (ST_MakeLine (geoma, geomb)::geography, geomp::geography), null from points;

-- PostGIS with the proposed algorithm 
/*  
     Run first the file "proposed_algorithm_postgis.sql" which is the implementation of the presented algorithm designed by 
	 Jose Martinez-Llario and Sergio Baselga.
	 It will add the stored procedure "STX_GeodesicMinDistance" that its needed in this code to get the closest distance from the geodesic line
*/  
insert into pointongeodesic (ncase, software, distance, geom) 
  select ncase,'PG_Paper', ST_Distance (geomp::geography, closestpoint::geography), closestpoint from 
    ( select STX_Geodesicmindistance (geoma, geomb, geomp) as closestpoint, ncase, geomp from points ) as foo;


--Insert distances from Oracle. 
--The results from ORacle are taken from the file geodesicmindistance_with_oracle.sql
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 1', 'Oracle', 23723.6561485326, null);
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 2', 'Oracle', 1013727.29054293, null);
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 3', 'Oracle', 3923508.12609771, null);

--Insert distances from Google Earth Engine. 
--The results from Google  are taken from the file geodesicmindistance_with_google_earthengine.js
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 1', 'Google', 23768.180843394905, null);
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 2', 'Google', 1015615.4291746469, null);
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 3', 'Google', 3929995.4526592367, null);

--Insert distances from ArcGIS (Near tool for geoprocessing after creating the geodesic line). 
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 1', 'ArcGIS', 23767.725756, 'SRID=4326;POINT (5.260427997 51.846089503)');
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 2', 'ArcGIS', 1010586.009156, 'SRID=4326;POINT (-21.93725208 54.928531473)');
insert into pointongeodesic (ncase, software, distance, geom) values ('Case 3', 'ArcGIS', 3928422.743555,'SRID=4326;POINT (18.349021882 37.978099032)');


	/* -- S.12345 (0.3mm in equador) S.123456789 (0.03um in equador) */
select ncase, software, distance, st_astext(geom) as degrees, (ST_AsLatLonText(geom, 'D°M''S.SSSSSSSSS"')) as DMS from pointongeodesic order by software, ncase;;

/*

	"ncase","software","distance","degrees","dms"
	"Case 1","ArcGIS","23767.725756","POINT(5.260427997 51.846089503)","51°50'45.922210800"" 5°15'37.540789200"""
	"Case 2","ArcGIS","1010586.009156","POINT(-21.93725208 54.928531473)","54°55'42.713302800"" -21°56'14.107488000"""
	"Case 3","ArcGIS","3928422.743555","POINT(18.349021882 37.978099032)","37°58'41.156515200"" 18°20'56.478775200"""
	"Case 1","Google","23768.1808433949","",""
	"Case 2","Google","1015615.42917465","",""
	"Case 3","Google","3929995.45265924","",""
	"Case 1","Oracle","23723.6561485326","",""
	"Case 2","Oracle","1013727.29054293","",""
	"Case 3","Oracle","3923508.12609771","",""
	"Case 1","PG_Paper","23767.7241838","POINT(5.26042849496107 51.8460892270512)","51°50'45.921217384"" 5°15'37.542581860"""
	"Case 2","PG_Paper","1010585.99883681","POINT(-21.9372910660488 54.9285314971169)","54°55'42.713389621"" -21°56'14.247837776"""
	"Case 3","PG_Paper","3928422.73531622","POINT(18.3490633150769 37.9781176779955)","37°58'41.223640784"" 18°20'56.627934277"""
	"Case 1","PostGIS","23768.1280775","",""
	"Case 2","PostGIS","1015610.50432715","",""
	"Case 3","PostGIS","3929987.74978078","",""


*/
	
select p0.ncase, p0.software, p1.software, p0.distance, p0.distance-p1.distance  from pointongeodesic p0, pointongeodesic p1 where p0.ncase = p1.ncase and p1.software = 'PG_Paper' order by p0.software, p0.ncase;
/*

	"ncase","software","software-2","distance","?column?"
	"Case 1","ArcGIS","PG_Paper","23767.725756","0.00157219999891822"
	"Case 2","ArcGIS","PG_Paper","1010586.009156","0.0103191899834201"
	"Case 3","ArcGIS","PG_Paper","3928422.743555","0.00823877984657884"
	"Case 1","Google","PG_Paper","23768.1808433949","0.456659594903613"
	"Case 2","Google","PG_Paper","1015615.42917465","5029.43033783685"
	"Case 3","Google","PG_Paper","3929995.45265924","1572.71734301653"
	"Case 1","Oracle","PG_Paper","23723.6561485326","-44.0680352674026"
	"Case 2","Oracle","PG_Paper","1013727.29054293","3141.29170612001"
	"Case 3","Oracle","PG_Paper","3923508.12609771","-4914.60921851033"
*/

/* Paper example

WITH minDistanceTable (x, p, ncase) AS (
  select STX_Geodesicmindistance (geoma, geomb, geomp), geomp, ncase from points 
)
SELECT ncase, ST_Distance (p::geography, x::geography) as mindist, 
              ST_AsLatLonText(x, 'D°M''S.SSSSSSSSS"') as DMS 
FROM minDistanceTable order by ncase;

"ncase","mindist","dms"
"Case 1","23767.7241838","51°50'45.921217384"" 5°15'37.542581860"""
"Case 2","1010585.99883681","54°55'42.713389621"" -21°56'14.247837776"""
"Case 3","3928422.73531622","37°58'41.223640784"" 18°20'56.627934277"""

*/



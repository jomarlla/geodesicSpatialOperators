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
  2 Run this file in the database.
  3 Run the file "GeodesicSpatialOp.sql" in the database to add the new spatial algorithms. 
  4 Run the file "geodesicintersection_with_new_paper_algorithm_postgis.sql" in the database.
  5 Run the file "test_A_B.sql" in the database to calculate the deviations from the tests
  
  The table "lines" contains the 8 geodetic lines to make the 4 intersection paper cases
  The table "inter" contains the PX1, PX2, PX3, PX4 coordinates calculated by Oracle, Google Earth, ArcGIS and PostGIS:
    
	If you want to get more information read these:
		Oracle Spatial: instructions in "geodesicintersection_with_oracle.sql"
		Google Earth Engine: instructions in "geodesicintersection_with_google_earthengine.js"
		ArcGIS and ArcGIS plus UTM: instructions in "geodesicintersection_with_arcgis.txt"
		PostGIS with local projections: lines 74-77 in this file
*/

create table lines (gid serial primary key, ncase varchar, name varchar, geom geometry (linestring, 4326));
insert into lines (ncase, name, geom) values ('Case 1', 'AB', 'SRID=4326;LINESTRING (14.5 54, 14.6 54.2)');
insert into lines (ncase, name, geom) values ('Case 1', 'CD', 'SRID=4326;LINESTRING (14.4 54.1, 14.7 54)');
insert into lines (ncase, name, geom) values ('Case 2', 'AB', 'SRID=4326;LINESTRING (5 52, 6 51.4)');
insert into lines (ncase, name, geom) values ('Case 2', 'CD', 'SRID=4326;LINESTRING (4.5 51.5, 5.5 52)');
insert into lines (ncase, name, geom) values ('Case 3', 'AB', 'SRID=4326;LINESTRING (29 42, -77 39)');
insert into lines (ncase, name, geom) values ('Case 3', 'CD', 'SRID=4326;LINESTRING (0 6, -22 64)');
insert into lines (ncase, name, geom) values ('Case 4', 'AB', 'SRID=4326;LINESTRING (-92 35, 52 40)');
insert into lines (ncase, name, geom) values ('Case 4', 'CD', 'SRID=4326;LINESTRING (20 -8, -95 49)');


--Add some special cases (reviewer comments)
insert into lines (ncase, name, geom) values ('Case 5', 'AB', 'SRID=4326;LINESTRING (5 90, 5 0)');
insert into lines (ncase, name, geom) values ('Case 5', 'CD', 'SRID=4326;LINESTRING (-30 70, 40 70)');
insert into lines (ncase, name, geom) values ('Case 6', 'AB', 'SRID=4326;LINESTRING (-175 80, 5 0)');
insert into lines (ncase, name, geom) values ('Case 6', 'CD', 'SRID=4326;LINESTRING (-30 60, 40 80)');
insert into lines (ncase, name, geom) values ('Case 7', 'AB', 'SRID=4326;LINESTRING (-170 85, 12 -15)');
insert into lines (ncase, name, geom) values ('Case 7', 'CD', 'SRID=4326;LINESTRING (-58 26, 120 75)');
insert into lines (ncase, name, geom) values ('Case 8', 'AB', 'SRID=4326;LINESTRING (105 63, 79 42)');
insert into lines (ncase, name, geom) values ('Case 8', 'CD', 'SRID=4326;LINESTRING (-167 38, -100 23)');
insert into lines (ncase, name, geom) values ('Case 9', 'AB', 'SRID=4326;LINESTRING (-42 40, 63 65.5)');
insert into lines (ncase, name, geom) values ('Case 9', 'CD', 'SRID=4326;LINESTRING (-41.8 40, 62.9 65.6)');



create table inter (gid serial primary key, ncase varchar, software name, geom geometry (point, 4326));


 
--The results for Oracle are taken from the file geodesicintersection_with_oracle.sql
--The results for Google  are taken from the file geodesicintersection_with_oracle_with_google_earthengine.js
--The results for ArcGIS/ArcGISUTM  are taken from the file geodesicintersection_with_arcgis.txt

insert into inter (ncase, software, geom) values ('Case 1', 'Oracle&Google', 'SRID=4326;POINT (14.528547835854896 54.05730103337946)');
insert into inter (ncase, software, geom) values ('Case 2', 'Oracle&Google', 'SRID=4326;POINT (5.22745698883988 51.8656617779676)');
insert into inter (ncase, software, geom) values ('Case 3', 'Oracle&Google', 'SRID=4326;POINT (-14.544131677975 54.6717159071832)');
insert into inter (ncase, software, geom) values ('Case 4', 'Oracle&Google', 'SRID=4326;POINT (-79.2501486314468 50.4442144186783)');
insert into inter (ncase, software, geom) values ('Case 1', 'ArcGIS', 'SRID=4326;POINT (14.528627872500067 54.057255744500083)');
insert into inter (ncase, software, geom) values ('Case 2', 'ArcGIS', 'SRID=4326;POINT (5.227448286000026 51.865640108000036)');
insert into inter (ncase, software, geom) values ('Case 3', 'ArcGIS', 'SRID=4326;POINT (-14.563892228999975 54.717012665500079)');
insert into inter (ncase, software, geom) values ('Case 4', 'ArcGIS', 'SRID=4326;POINT (-79.282819925499950 50.479069137000067)');
insert into inter (ncase, software, geom) values ('Case 1', 'ArcGISUTM', 'SRID=4326;POINT (14.5285483925224 54.0573011622007)');
insert into inter (ncase, software, geom) values ('Case 2', 'ArcGISUTM', 'SRID=4326;POINT (5.227448471388889 51.865640713611114)');
insert into inter (ncase, software, geom) values ('Case 3', 'ArcGISUTM', 'SRID=4326;POINT (-14.563855968611112 54.71702961777778)');
insert into inter (ncase, software, geom) values ('Case 4', 'ArcGISUTM', 'SRID=4326;POINT (-79.2828200318611 50.479069611944446)');

-- PX1,PX2,PX3,PX4 with PostGIS transformation to local projection
insert into inter (ncase, software, geom) select 'Case 1','PostGIS', ST_Intersection (l1.geom::geography, l2.geom::geography)::geometry from lines l1, lines l2 where l1.name = 'AB1' and l2.name = 'CD1';
insert into inter (ncase, software, geom) select 'Case 2','PostGIS', ST_Intersection (l1.geom::geography, l2.geom::geography)::geometry from lines l1, lines l2 where l1.name = 'AB2' and l2.name = 'CD2';
insert into inter (ncase, software, geom) select 'Case 3','PostGIS', ST_Intersection (l1.geom::geography, l2.geom::geography)::geometry from lines l1, lines l2 where l1.name = 'AB3' and l2.name = 'CD3';
insert into inter (ncase, software, geom) select 'Case 4','PostGIS', ST_Intersection (l1.geom::geography, l2.geom::geography)::geometry from lines l1, lines l2 where l1.name = 'AB4' and l2.name = 'CD4';




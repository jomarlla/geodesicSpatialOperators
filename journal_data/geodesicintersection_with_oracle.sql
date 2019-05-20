--Tested with:
--Oracle Database 12c Enterprise Edition Release 12.2.0.1.0 - 64bit Productionsystems
--Oracle Database 18c Enterprise Edition Release 18.0.0.0.0 - Production Version 18.3.0.0.0

--You must open a connection to any database and run the followings SQL commands:

CREATE TABLE LINES (
  NAME VARCHAR(3) PRIMARY KEY,
  GEOM SDO_GEOMETRY);
  
INSERT INTO USER_SDO_GEOM_METADATA (TABLE_NAME, COLUMN_NAME, DIMINFO, SRID) 
    VALUES ('LINES','GEOM', 
       MDSYS.SDO_DIM_ARRAY 
        (MDSYS.SDO_DIM_ELEMENT('Longitude', -180, 180, 0.05), --Tolerance = 0.05 (*)
         MDSYS.SDO_DIM_ELEMENT('Latitude', -90, 90, 0.05) ),  --Tolerance = 0.05 (*)
       NULL);
	   
-- (*) SDO_TOLERANCE NUMBER, 
-- If a function accepts an optional tolerance parameter and this parameter is null or not specified, the SDO_TOLERANCE value of the layer is used.
-- https://docs.oracle.com/cd/B14117_01/appdev.101/b10826/sdo_intro.htm#i884635
-- For Geodetic Data the tolerance is in meters. The minimum tolerance allowed by Oracle is 0.05.


--Case 2
INSERT INTO LINES VALUES(
  'AB2',
  SDO_GEOMETRY(
    2002,
    4326, -- EPSG:4326  = WGS84
    NULL,
    SDO_ELEM_INFO_ARRAY(1,2,1), -- LineString Type
    SDO_ORDINATE_ARRAY(5,52,6,51.4)
  )
);

INSERT INTO LINES VALUES('CD2', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(4.5,51.5,5.5,52) ) );

--Case 3
INSERT INTO LINES VALUES('AB3', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(29,42,-77,39) ) );
INSERT INTO LINES VALUES('CD3', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(0,6,-22,64) ) );

--Case 4
INSERT INTO LINES VALUES('AB4', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(-92,35,52,40) ) );
INSERT INTO LINES VALUES('CD4', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(20,-8,-95,49) ) );

--Case 1
INSERT INTO LINES VALUES('AB1', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(14.5,54,14.6,54.2) ) );
INSERT INTO LINES VALUES('CD1', SDO_GEOMETRY( 2002, 4326, NULL, SDO_ELEM_INFO_ARRAY(1,2,1), SDO_ORDINATE_ARRAY(14.4,54.1,14.7,54) ) );


SELECT sdo_util.to_wktgeometry ( sdo_geom.sdo_intersection (l1.geom, l2.geom, 0.05)) FROM lines l1, lines l2 WHERE l1.name = 'AB1' AND l2.name = 'CD1';
SELECT sdo_util.to_wktgeometry ( sdo_geom.sdo_intersection (l1.geom, l2.geom, 0.05)) FROM lines l1, lines l2 WHERE l1.name = 'AB2' AND l2.name = 'CD2';
SELECT sdo_util.to_wktgeometry ( sdo_geom.sdo_intersection (l1.geom, l2.geom, 0.05)) FROM lines l1, lines l2 WHERE l1.name = 'AB3' AND l2.name = 'CD3';
SELECT sdo_util.to_wktgeometry ( sdo_geom.sdo_intersection (l1.geom, l2.geom, 0.05)) FROM lines l1, lines l2 WHERE l1.name = 'AB4' AND l2.name = 'CD4';

--Case 1: POINT (14.5285478358547 54.0573010333798) ->      14° 31' 42.77221    54° 3' 26.28372
--Case 2: POINT (5.22745698883988 51.8656617779676) -> 		5° 13' 38.84516 	51° 51' 56.38240
--Case 3: POINT (-14.544131677975 54.6717159071832) -> 		-14° 32' 38.87404	54° 40' 18.17727
--Case 4: POINT (-79.2501486314468 50.4442144186783) -> 	-79° 15' 0.53507	50° 26' 39.17191
--Comparison with Google Earth Engine: (1E-12 difference) around 1 micrometer on the ecuador
		--Case 1: 14.528547835854896,54.05730103337946
		--Case 2: 5.227456988839927,51.865661777967624
		--Case 3: 14.544131677975026,54.67171590718321
		--Case 4: 79.25014863144679,50.44421441867829

--sdo_intersection - tolerance, for geodesic data is in meters and 0.05 is the minimum amount in Oracle.
--https://docs.oracle.com/database/121/SPATL/data-model.htm#SPATL450


-- drop table lines;
-- delete from user_sdo_geom_metadata where table_name = 'LINES';
 
select sdo_geom.sdo_length(geom) from lines;



// Google earth enginee
// Last version tested from 

// Run this code into the google online editor:
// https://code.earthengine.google.com/#

//WGS84 by default

/*  https://developers.google.com/earth-engine/api_docs#

	ee.Geometry.Point.distance (right, maxError, proj)
	
	maxError	ErrorMargin, default: null	
	The maximum amount of error tolerated when performing any necessary reprojection.

	proj	Projection, default: null	

	The projection in which to perform the operation. If not specified, the operation will be performed in a 
	spherical coordinate system, and linear distances will be in meters on the sphere.
*/

//Checking point to line minimum distance (the two first cases are AB2 and AB3 from the previous example).
var AB1 = ee.Geometry.LineString ( [[5, 52], [6, 51.4]] );
var AB2 = ee.Geometry.LineString ( [[29, 42], [-77,39]] );
var AB3 = ee.Geometry.LineString ( [[29, 42], [-70, -35]] );

var X1 = ee.Geometry.Point ( [5.5, 52] );
var X2 = ee.Geometry.Point ( [-22, 64] );
var X3 = ee.Geometry.Point ( [-22, 64] );

print (AB1.distance (X1, 0.001));
print (AB2.distance (X2, 0.001));
print (AB3.distance (X3, 0.001));

//Case 1: 23768.180843394905
//Case 2: 1015615.4291746469
//Case 3: 3929995.4526592367
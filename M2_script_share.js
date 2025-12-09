//////////////////////////////////////////import NAIP and filter by sites
//Create a variable called NAIP to point to the NAIP data collection
var NAIP = ee.ImageCollection('USDA/NAIP/DOQQ');

var bound = table
  .filterBounds(geometry2); //filter with defined filter location

Map.addLayer({
  eeObject:bound,
  shown:false
});

//Filter NAIP imagery and assess start and end dates
var naipImage1 = NAIP
  .filterBounds(bound) //filter with defined filter location
  .filterDate('2006-01-01','2020-12-30') //Limit images to the ones collected in 2010
  .select(['R', 'G', 'B'])
  .sort("system:time_start");
  //.mosaic(); //Convert selected images into one composite image rather than an image collection.
var i12 = naipImage1.map(function(img){return img.clip(bound)});
print(i12);

var Dates = i12.aggregate_array('system:time_start');
//print('Dates Acquired:', Dates);
print(Dates.map(function(da){return ee.Date(da)}));

var Dates = i12.aggregate_array('system:time_end');
//print('Dates Acquired:', Dates);
print(Dates.map(function(da){return ee.Date(da)}));
////////////////////////////////////////


/////////////////////////////////function below named "proce" is for running NAIP image clipping
/////////////////////////////////entropy calculation, Geary's C, and final tree vs. shrub classification.
var proce = function(x,y,z,name){
    //Filter NAIP imagery
  var naipImage = NAIP
    .filterBounds(bound) //filter with defined filter location
    .filterDate(x,y) //Limit images to the ones collected in 2010
    .sort("system:time_start")
    .mosaic(); //Convert selected images into one composite image rather than an image collection.
  
  Map.addLayer({
    eeObject:naipImage.clip(bound),
    shown:true
  });
  
  Export.image.toDrive({
    image: naipImage,
    region: bound,
    description: name,
    fileNamePrefix: name,
    scale: 1
  });
  
  Map.addLayer({
    eeObject:table,
    shown:false
  });
  Map.centerObject(bound);
  //print(table);
  
  
  ///////////////// compute entropy (texture)
  // Load a high-resolution NAIP image.
  //var image = ee.Image('USDA/NAIP/DOQQ/m_4207710_sw_18_1_20110603');
  var image = naipImage.clip(bound);
  
  // Get the red band.
  var nir = image.select('R');
  
  // Define a neighborhood with a kernel.
  var square = ee.Kernel.square({radius: 4});
  
  // Compute entropy
  var entropy = nir.entropy(square);
  
  // Create a list of weights for a 9x9 kernel.
  var list = [1, 1, 1, 1,1];
  
  // The center of the kernel is zero.
  var centerList = [1, 1, 0, 1, 1];
  
  // Assemble a list of lists: the 9x9 kernel weights as a 2-D matrix.
  var lists = [list, list, centerList, list, list, ];
  
  // Create the kernel from the weights.
  // Non-zero weights represent the spatial neighborhood.
  var kernel = ee.Kernel.fixed(5, 5, lists, -4, -4, false);
  
  // Convert the neighborhood into multiple bands.
  var neighs = nir.neighborhoodToBands(kernel);
  
  //////////////// Compute local Geary's C (texture) and use it to filter for connected pixels
  // Geary's C is used to help improve Tree classification based on entropy results
  var gearys = nir.subtract(neighs).pow(2).reduce(ee.Reducer.sum())
              .divide(Math.pow(5, 2));
  
  //selected connected pixels that are greater than a certain number of pixels (patchsize)
  //arrgresive on the first filter, then narrow down with second filter
  //otherwise, end results will not have connected pixels
  
  //first compute filter out local Geary's C greater than z value 
  //count number of connected pixels smaller than 100
  var patchsize = gearys.lte(z).connectedPixelCount(100, false);
  var gearys_processing = gearys.lte(z).multiply(patchsize.gte(100)).eq(0);
  
  //second compute and filter out local Geary's C greater than z value 
  //count number of connected pixels smaller than 30
  var patchsize1 = gearys_processing.eq(1).connectedPixelCount(30, false);
  
  //filter out patchsize with less than 30 as trees
  //set tree pixel values to 1, and set non-tree pixel to 0
  var gearys_processing1 = gearys_processing.eq(1).multiply(patchsize1.gte(30)).eq(0);
  var gearys_processing2 = gearys_processing1.add(1);
  
  //Export final classified raster image to Google Drive
  Export.image.toDrive({
    image: gearys_processing2,
    region: bound,
    description: name,
    fileNamePrefix: name,
    scale: 1
  });
}

/////////////////////////////////////////////////////code below run the entire "proce" function
/////////////////////////////////////////////////////function inputs are start date, end date, 
/////////////////////////////////////////////////////entropy threhold, and final output image "label"
proce('2017-01-01','2017-12-31',200,"2017")

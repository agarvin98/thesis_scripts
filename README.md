# thesis_scripts
Drivers of Community Assembly in Recently Disturbed Habitat -- all code.

Helllllooooo!!!!!

The scripts in this sub-folder labelled 1-6 show the iterative steps towards building a species distribution model. Always double check model output at each step, especially if you are running a loop to check if it works properly first. 

1. Raw zerofill -- filtering data
2. Swiftlet thin plot -- spatial autocorrelation fixing biases
3. Layer processing -- cropping, buffering, projecting, masking layers..
4. Select07GLM -- Removing correlated variables 
5. ENMEval NullMods -- ENMEval, model selection, model evaluation, null model iterations, niche breadth and overlap calculations. 
6. Sabah processing script -- Selection of polygons based off suitability across each species and accessibility. For field research and data collection with Survey123 and ArcFieldMaps. Links below:

Links for surveys and details:
Draft of some public access surveys I made for (a) swiftlet presence information, so taking gps points and collecting data with embedded photos for species presence/nest/foraging. Users (us) will repeat the survey for each species/colony found, it will associate it with the GPS point where you collected the data. (b) will ideally be done on a device that is open for the duration of our sampling where we can also track the survey path around our buffers. [subject to change]

      Swiftlet presence/nests/foraging survey: https://arcg.is/1Hq8nv0
      Swiftlet environmental characteristics survey: https://arcg.is/0vaTPP
      This is a draft attempt for accessing and locating centroids of polygons with 200m buffer:  https://fieldmaps.arcgis.app/itemID=c981a65078d943b5a1241ec9ab649e2f&referenceContext=open&portalURL=https%3A%2F%2Fpomona.maps.arcgis.com
      

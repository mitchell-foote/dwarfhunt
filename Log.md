# Log for my studies

## Sprint Week

- In order to get a random set of data in a range, use the np.random.uniform, which has the values for min and max, instead of random.random, which doesn't

- In order to do np arrays like js maps, best bet is to try a lambda function, and vectorize it beforehand if you need to. You don't always need to. I do wish they had a map function.

- Floating point number issues are usually why you see larger clifs in the graph, not a data problem. 

- Make sure to import all of the sub stuff you need for species package

- Brown dwarfs (especally in the L/T transition) don't always show cooler = redder. Around T I discovered that methane can form in large enough quantities to block some of the H filter I was using from the tutorial. 
  
- Reddening vector shows the direction that should occur due to to adding dust to the measured values.

- If you get a dwarf that isn't matching the reddening factor you're tracking, it's probably not that factor that's causing the issue. 

- Best place to check the filter names is https://svo2.cab.inta-csic.es/theory/fps/index.php Don't expect them to just be their name.

- There's a difference between the filter value of the filter only, and the whole system. The whole system contains information about QE, throughput, as well as the filter values. 

-  Filter values come back as a 2d array of wavelengths and transmissions

- Working through repro-ing, it's time to start learning about atmo, and how it impacts brown dwarfs and the measurements. 

- Adding new variables in models increases the time to calclate the various grid points, so adding a new variable instead of using a constant requires a significant scientific reason to do so. 

- Another key point is that there isn't a magic function that can be called to calculate the perfect values in a downloadable model. To get around this scientists check certain points, and then have functions interpolate around them. Is there a function that could make that work? Yes. But running that function is for super-computers, not for python. 

- K_zz stands for vertical atmospheric mixing. A high mixing value means the atmo has a larger varying chem, and can produce unexpected spectrum. 
  
- Download a full model causes a huge memory usage. Not really what you want to do, so you need to set a wavel_range and resolution, so that it cuts the amount of data that you need to work with. For what I'm doing, I don't need the whole thing. 

- I've noticed from clauding around that the elf-owl models don't go as high as the JWST filter. I need to go through and figure out how the paper dealt with that kind of infomration. 

- A good tutorial for finding synthetic photometry lives here: https://species.readthedocs.io/en/latest/tutorials/model_spectra.html It's not typesafe, but it works

- To get flux, you want to use the ReadModel.get_flux functions to pull that out. Specify which filter value so you can do a full color - color fit. 

- To get Mag, you want to use the same class but get_magnitude

- Sigh... time to remember that abs mag and apparent mag are what would be if the object was viewed from 10 pc away, which turns out to be the same if you specify 10 pc. Huzzah, not crazy for seeing the same number twice. 



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

- We're running into an issue where the model specified in the paper don't go through the final F1550C filter's effective wavelength. elf-owl ends at 15um, while the F1550C's effective range is (np.float64(14.912), np.float64(16.219)), so barely hitting that lower end, but not getting anywhere close the peak of that filter. I get an error on species when attempt to load the model for those wavelengths. 

- This might be a good example where they're doing some sort of custom model entry, or I'm missing a species param to force it to extrapolate, but I'm not seeing it on my end. I've also checked the original elf-owl paper (https://arxiv.org/pdf/2402.00756) to verify that the wavelength values reported by species are correct. From what I can tell, they are. 

- But looking at the paper Michelson et al 2025, it says they used the library to get that info. Looks like a good time to reach out to Zac / the paper authors to see how they got those values. 


- In order to unblock myself, I'm going to try and create the graph with simulated photometry using the sonora-diamondback-highres model. It's smaller, and doesn't have the same, or as many params, so you won't be able to get the same resolution, but at least we'll have enough for the color/color diagram


- After doing some looking, I was able to determine that the diamondback model isn't a good fit because it doesn't focus on the cold end of the spectrum. It's not really what I'm looking for here, because it doesn't fit the science I'm testing. 


## 2026-07-27

### Staring the week
- We've got the graph as far as I can go without additional information. 
- I've reached out to Ashley Messier with a request to get the galaxy information and how to understand how they got the final JWST filter
- It's time to start jumping on some of the technical learnings that I'll need as well as do some galaxy exploration. Maybe I can find the data in other places. 

### What's blocked, and on whom
- Blocked on getting the filter data, as well as the Galaxy data. 

### What I learned 
- One of the issues that I ran into while working on the graph is that we don't have data for some of the points. So we have a grid, but not every point on the grid is filled. 
- To combat that, we look at the data and only sample from the points we know we have. It makes the graph look segmented, but we can play around with it after I get some of the other data. 

### Next
- Start on the stats + astro work that I generated to help get me up to speed. 

## 2026-07-29

Heard back from Ashley! She was very helpful and gave some advice on how to get the galaxy models, as well as some pointers for where my exploration can go later in the fall. 

Basically, she's going to reach out to the PI to see how they were able to get the F1550C's filter data to work with the elf-owl models, but also gave me some excellent ideas on how to extend my work by verifying the spread factor on a lot of the new models. Facinating stuff, but stuff that I'll need to check out after I do some replication. 

The specific model she's referening is the brand new exo-remk26 model that just came out this year. The data can be found here: https://lesia.obspm.fr/exorem/YGP_grids/Exo-REMk26/High_Res_grid_2026/R200k_cloudyfsed_2026/

I've also been working on making the galaxy graph, and I've realized that it's often helpful to look at the original author's website instead of a mirror, as the mirror doesn't have to have the complete data set. I ran into this with the SWIRE subset, and the original author's (Polletta) distro. 

I did some research into how Vega can affect the magnitude readings that come back from species. Thankfully this filters out when doing a color-color analysis, as everything is being affected by the same Vega drops in MIR range. However, when you redshift too far, you might run into issues with magnitude adjustments. 

Also dug into the difference between F_λ and νF_ν, and why astronomers like the log log adjusted graphs in νF_ν, as it shows it adjusted for the log scale. 



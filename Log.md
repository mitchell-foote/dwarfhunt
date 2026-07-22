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


  

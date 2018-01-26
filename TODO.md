* Start up with the Sun and Jupiter in the correct places and get them moving correctly
  - Ignore the impacts the asteroids can have on their orbit, probably
  - Think about getting different timesteps for different asteroids...
  
* Change the initialization for the asteroids to have a certain semi-major axis & eccentricity and then pick a random starting point along the orbit (with matching speeds, obv)
  - Initially have no eccentricity
  - Probably want to just use x&y location&speeds and use a small random range for the z-speeds, although that is further down the line
    
* Find an easy way to throw out asteroids which become too elliptical or go too far out (need2track changes in eccentricity for the first which might be not-practical...)

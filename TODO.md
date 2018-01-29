* ~~Expand on the Kirkwoods class so each of them is aware *where* the sun and jupiter are, rather than having to figure this out each time again for each asteroid and pass it as an argument~~
  - ~~Once this is implemented the threebody simulation can easily become the nbody sim, as it simply becomes updateplanet into updateasteroids~~

* Change the initialization for the asteroids to have a certain semi-major axis & eccentricity and then pick a random starting point along the orbit (with matching speeds, obv)
  - Current working version puts them at the x=0 and y=max location (with matching speed)
  - Initially have no eccentricity
  - Probably want to just use x&y location&speeds and use a small random range for the z-speeds, although that is further down the line
    
* Find an easy way to throw out asteroids which become too elliptical or go too far out (need2track changes in eccentricity for the first which might be not-practical...)
  - As asteroids are now aware how far from the sun they are, this should be an easy thing to implement
  - Might be an idea to also throw away asteroids which become too attached to Jupiter (moon-like)

* Find a way to have the asteroids be aware what their period is
  - Possibly track how long it takes them to cross from negative x -> positive x? Assuming they are following (mostly) circular orbits that should work
  - Can even track eccentricity by checking the difference in positive x time versus negative x time which might allow for throwing out very eccentric orbits if we want to

* Later additional options:
  - ~~Think about getting different timesteps for different asteroids...~~ It's probably easier to ask all asteroids how small it wants the timestep to be and use the smallest timestep for all objects, rather than have different timesteps for different objects at a certain time

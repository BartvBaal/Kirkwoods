* ~~Expand on the Kirkwoods class so each of them is aware *where* the sun and jupiter are, rather than having to figure this out each time again for each asteroid and pass it as an argument~~
  - ~~Once this is implemented the threebody simulation can easily become the nbody sim, as it simply becomes updateplanet into updateasteroids~~

* Change the initialization for the asteroids to have a certain semi-major axis & eccentricity and then pick a random starting point along the orbit (with matching speeds, obv)
  - ~~Current working version puts them at the x=0 and y=max location (with matching speed)~~
  - ~~Initially have no eccentricity~~ Eccentricity is an option for asteroids, atm always 0
  - ~~Probably want to just use x&y location&speeds and use a small random range for the z-speeds, although that is further down the line~~ 3D transition done, ~~currently no z speeds~~
    
* Find an easy way to throw out asteroids which become too elliptical or go too far out (need2track changes in eccentricity for the first which might be not-practical...). ~~Change in setup no longer throws away far-off asteroids, need to re-establish a way to do this (WIP).~~
  - ~~As asteroids are now aware how far from the sun they are, this should be an easy thing to implement.~~ Currently throw away all asteroids which go further than ~~6.5AU~~ 7 AU. Considering throwing away close-to-sun asteroids (within 1.2 AU??) too.
  - Might be an idea to also throw away asteroids which become too attached to Jupiter (moon-like)

* ~~Find a way to have the asteroids be aware what their period is~~ First version live; currently assumes P² = a³ ~~and checks if the position one full period back is within 10% distance of the final position. Need to figure out to track backwards ?~~
  - ~~Possibly track how long it takes them to cross from negative x -> positive x? Assuming they are following (mostly) circular orbits that should work~~
  - ~~Can even track eccentricity by checking the difference in positive x time versus negative x time (and the same for y and z possibly) which might allow for throwing out very eccentric orbits if we want to~~
  - Got reminded by astrophysicists that there is a formula to calculate this. Adapted final plots for this too, need to discuss using this method to throw out runaway asteroids.

* ~~Initial time-check shows get_distance takes up a lot of time - would be nice if this can be optimized. Most other time goes into update_asteroid, also a point for optimization~~
  - ~~Going to attempt to speed this up by using numpy/scipy for calculations~~ Praise be numpy array
  - ~~Can solve part of the get_distance issue by just putting the Sun at the center and ignoring the whole center-off-mass @ center part (pretend the sun is the center of mass). This would also slightly cut down on the runtime of update_planet~~ Currently not an issue

* Later additional options:
  - ~~Think about getting different timesteps for different asteroids...~~ It's probably easier to ask all asteroids how small it wants the timestep to be and use the smallest timestep for all objects, rather than have different timesteps for different objects at a certain time

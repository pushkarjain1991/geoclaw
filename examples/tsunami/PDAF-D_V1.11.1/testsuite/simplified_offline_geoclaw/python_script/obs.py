import numpy as np

def make_obs(dimx, dimy, dxobs, dyobs,stddev_obs,truefield, testing = False):

    if testing:
        obs_error = 0.0
    else:
        obs_error = stddev_obs*np.random.randn(dimx,dimy)
    
    full_obs = truefield + obs_error
    obs = np.zeros_like(truefield)- 999
    obs[dxobs:dimx:dxobs,dyobs:dimy:dyobs]=full_obs[dxobs:dimx:dxobs,dyobs:dimy:dyobs]

    print obs
    np.savetxt("obs.txt",obs)

if __name__=='__main__':
    make_obs(np.zeros((51,51)))

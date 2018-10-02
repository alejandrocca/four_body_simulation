def expand_body(newsys, oldcom, syst):

    for mem in newsys:
        mem.position += oldcom.position
        mem.velocity += oldcom.velocity
       
    ## removing fictional particle
    ## first check if m1+m2 = p.mas
    if(newsys.mass.sum() != oldcom.mass):
        raise ValueError("components mass not equal to center of mass\n{:s}!={:s}".format(newsys.mass.sum(), oldcom.mass))
    syst.remove_particle(oldcom)
    syst.add_particles(newsys)
   
    return sys

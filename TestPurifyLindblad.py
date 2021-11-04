import numpy as np
import pyten as ptn

#Parameters
N_sites=2
max_bosons=5

#Lattice
ferm_bos_sites = []
for i in range(N_sites):
    ferm_bos_sites.append(1) #fermionic sites
    ferm_bos_sites.append(1) #purified fermionic sites
    ferm_bos_sites.append(0) #bosonic sites
    ferm_bos_sites.append(0) #purified bosonic sites

lat = ptn.mp.lat.u1u1u1.genFermiBoseLattice(ferm_bos_sites,max_bosons) #,purify=True  #length: 2*N_sites
##lat.save("lat_test")
pp_lat = ptn.mp.proj_pur.proj_purification(lat, [0], ["a", "ah"])


#Build infiite temperature state
def GenMaxEntStateBosons(site,pp_lat,max_bosons):
    c = [pp_lat.get('I')]
    d = pp_lat.get('I')
    for i in range(max_bosons):
        d *=1./(i+1)*pp_lat.get('a',site+1)*pp_lat.get('a',site+3)*pp_lat.get('ah',site)*pp_lat.get('ah',site+2)
        c += [d]
    c = ptn.mp.addLog(c)
    return c

def GenMaxEntStateFermions(site,pp_lat,max_bosons):
    c = [pp_lat.get('I')]
    c += pp_lat.get('chu',site)*pp_lat.get('chu',site+1)
    c += pp_lat.get('chd',site)*pp_lat.get('chd',site+1)
    c += pp_lat.get('chu',site)*pp_lat.get('chu',site+1)*pp_lat.get('chd',site)*pp_lat.get('chd',site+1)

    c = ptn.mp.addLog(c)
    return c


#Generate maximally entangled state from the vacuum
Op = pp_lat.get('I')
for i in range(N_sites):
    print(i,)
    #Op *= GenBathThermalState(1+i*(4*N_bos_sites+1)+4*j,pp_lat,max_bosons,beta,bos_freq_vec[j])
    #Op.truncate()

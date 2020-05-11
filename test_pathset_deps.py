'''
Test suite for functions I have written that the Pathset module depends on
'''
import pytest
from calc_aoi import slw2aoi, get_rayparam
from sphe_trig import vincenty_dist, vincenty_direct
from numpy import around
from obspy.clients.iris import Client

@pytest.mark.parametrize("d, slw, expect",
                         [(0, 0, 0), (50,200,8.16098368022881)])
def test_slw2aoi(d,slw,expect):
    ''' Tests the slw2aoi function '''
    output = slw2aoi(d,slw)
    
    assert output == expect

@pytest.mark.parametrize("evdp, dist, ph, expect",
                         [(10,100,'SKS',281.70596), (10,100,'SKKS',413.07653), 
                         (300,75,'ScS',461.90343)])
def test_get_rayparam(evdp, dist, ph, expect):
    out = get_rayparam(evdp, dist, ph)
    output =  around(out, decimals=5)
    assert output == expect
@pytest.mark.parametrize("lat1,lon1,lat2,lon2",
                        [(5,5,10,10), (25,12,33,45), (0,0,12,-35),(-70,-120,-45,-20)])
def test_vincenty_dist(lat1,lon1,lat2,lon2):
    testClient = Client()
    distaz = testClient.distaz(lat1,lon1,lat2,lon2)
    C_d,C_baz =  distaz['distance'],distaz['backazimuth']
    d,a = around(vincenty_dist(lat1,lon1,lat2,lon2),decimals=4)
    dd = abs(C_d - d)
    da = abs(C_baz - a)
   
    assert [dd, da] <= [1e5, 1e5]
    
@pytest.mark.parametrize("lat1,lon1,azi,dist,expect",
                        [(5, 5, 44.64334, 7.03999,[10.0, 10.0])])
def test_vincenty_direct(lat1,lon1,azi,dist,expect):
    lat2,lon2 = vincenty_direct(lat1,lon1,azi,dist)
    output = [lat2,lon2]
    
    assert output == expect 


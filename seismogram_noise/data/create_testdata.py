import instaseis
db = instaseis.open_db('syngine://ak135f_2s')
src = instaseis.Source(m_rr=1e10, latitude=0.0, longitude=0.0, depth_in_m=10e3)
rec = instaseis.Receiver(latitude=50, longitude=8, network='YY', station='TEST')
st = db.get_seismograms(src, rec, dt=0.1, kind='velocity')
print(st)
st.write('testdata.mseed', format='MSEED')

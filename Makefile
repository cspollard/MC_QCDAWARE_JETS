RivetMC_QCDAWARE_JETS.so: MC_QCDAWARE_JETS.cc
	rivet-buildplugin RivetMC_QCDAWARE_JETS.so MC_QCDAWARE_JETS.cc `fastjet-config --prefix`/lib/libQCDAware.a

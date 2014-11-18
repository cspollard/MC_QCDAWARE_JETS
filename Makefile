RivetMC_QCDAWARE_JETS.so: MC_QCDAWARE_JETS.cc FinalPartons.cc FinalPartons.hh
	rivet-buildplugin RivetMC_QCDAWARE_JETS.so MC_QCDAWARE_JETS.cc FinalPartons.cc `fastjet-config --prefix`/lib/libQCDAware.a

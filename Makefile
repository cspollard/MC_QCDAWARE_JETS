RivetMC_QCDAWARE_JETS.so: MC_QCDAWARE_JETS.cc UserInfoParticle.hh
	rivet-buildplugin RivetMC_QCDAWARE_JETS.so MC_QCDAWARE_JETS.cc `fastjet-config --prefix`/lib/libQCDAwarePlugin.a

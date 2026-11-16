# KVASIR
Proposal for high resolution neutron spectrometer for the European Spallation Source

## Instriment concept
KVASIR is an indirect time of flight instrument concept designed for high resolution neutron spectroscopy on single crystals. The instrument aims at an elastic energy resolution of ~4$\mu$eV FWHM and a momentum transfer resolution of 0.03Å$^{-1}$. 

The instrument offers two modes of operation. The "high resolution" mode is the one described above. The second mode is using second order Bragg scattering of the analyzer. This mode provides an enhanced dynamic range of the instrument and will therefore be refered to as the "wide dynamic range" mode. 

### Lore
Kvasir is the norse god of wisdom. He was killed by evil dwarfs how brewed the first mjod on his blood. The point being, that drinking mjod will make you infinitely wise.

## Running the code 
Instrument is based on the McStas package.
Virtual experiements are run using master.py with options:
- experiment: incoherent, powder, magnon, single_crystal
- mode: high_res, wide_dyn_range

Called like this: python3 master.py incoherent high_res

Simulations of incoherent and powder experiments takes around 4 hour to complete using my Mac(Silicon), The magnon and single_crystal experiments are included, but running them locally is not recomended.

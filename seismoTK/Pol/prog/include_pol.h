cccc 
cc    parameter (maxns=65536,maxtr=1,maxfr=121)
c TFR:  use nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
cc    parameter (nstf=65536,nftf=32769)
c nax>= nsmp >= maxns
cc    parameter (nax=65536)

cccc 32768  works fine on my new laptop
cc      parameter (maxns=32768,maxtr=1,maxfr=321)
c TFR:  use nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
cc      parameter (nstf=32768,nftf=16385)
c nax>= nsmp >= maxns
cc      parameter (nax=32768)

cccc 16384
cc      parameter (maxns=16384,maxtr=1,maxfr=101)
c TFR:  use nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
cc      parameter (nstf=16384,nftf=8193)
c nax>= nsmp >= maxns
cc      parameter (nax=16384)

cccc 8192 (SAFE) works fine on my centrino laptop with 1Gb memory! 
      parameter (maxns=8192,maxtr=1,maxfr=1024)
c TFR:  use nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
      parameter (nstf=8192,nftf=4097)
c nax>= nsmp >= maxns
      parameter (nax=8192)

cccc 4096
cc     parameter (maxns=4096,maxtr=1,maxfr=1021)
c TFR:  use nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
cc     parameter (nstf=4096,nftf=2049)
c nax>= nsmp >= maxns
cc     parameter (nax=4096)


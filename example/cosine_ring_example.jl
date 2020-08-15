using VIDA

cd(@__DIR__)
img_el1 = load_ehtimfits("elliptical_gaussian_rot-45.00m87Scale_seed_23_simobs_netcal_scanavg-z0.6-s100-t0-v0-l0-p50-e0.000.fits")
plot(img_el1)
bh_el1 = Bhattacharyya(img_el1)

img_el2 = load_ehtimfits("elliptical_gaussian_rot-45.00m87Scale_seed_23_simobs_netcal_scanavg-z0.6-s100-t1-v100-l10-p40-e0.000.fits")
plot(img_el2)
bh_el2 = Bhattacharyya(img_el2)

#GRMHD imgs
img_gr2 = load_ehtimfits("/home/ptiede/Research/Projects/Waffle/WaffleTest/Images/3597/dir_001/Good/image_a-0.94_1242_17_0_230.e9_6.2e9_8.03224e+28_20_SANE_3597_scale-1.1_seed_14_simobs_netcal_scanavg-z0.8-s100-t0-v1-l10-p50-e0.020.fits")
plot(img_gr2)
bh_gr2 = Bhattacharyya(img_gr2)

img_gr1 = load_ehtimfits("/home/ptiede/Research/Projects/Waffle/WaffleTest/Images/3597/dir_001/Good/image_a+0.5_0007_163_0_230.e9_6.2e9_1.59542e+25_40_MAD_3597_scale-1.0_seed_21_simobs_netcal_scanavg-z0.6-s100-t0-v0-l100-p50-e0.000.fits")
plot(img_gr1)
bh_gr1 = Bhattacharyya(img_gr1)
#Elliptical ring example runs

#Start with a general gaussian ring filter
θg = GeneralGaussianRing(rand(8)) + 1.0*Constant()
lowerg = [1.0  ,0.1, 0.001, 0.0, 0.001, -π, -60.0, -60.0, 1e-6]
upperg = [30.0,20.0, 0.999,   π, 0.999,  π,  60.0,  60.0, 1]

resg_e1 = bbextract(bh_el1,
                    θg,
                    lowerg,
                    upperg;
                    TraceMode=:compact,
                    MaxFuncEvals=20000)
p=triptic(img_el1, resg_e1[1])
display(p)

resg_e2 = bbextract(bh_el2,
                    θg,
                    lowerg,
                    upperg;
                    TraceMode=:compact,
                    MaxFuncEvals=20000)
p = triptic(img_el2, resg_e2[1])
display(p)

#Okay that's cool and all but let's say I actually want to model
#the weird blobby structure of the image well now we can do that
θ23 = CosineRing{2,3}(20.0,
                      [5.0,1.0],
                      [0.0],
                      0.2, 0.0,
                      [0.5, -0.5, 0.1],
                      [π/6, 0.0, 0.0],
                      -10.0, 10.0
                      ) + 1.0*Constant()
lower23 = [1.0,
           0.1, -5.0,
           -π,
           0.01, 0.0,
           0.01,-0.99,-0.99,
           -π,-π,-π,
           -60.0, -60.0,
           1e-6]

upper23 = [40.0,
           15.0, 5.0,
           π,
           0.99,π,
           0.99,0.99,0.99,
           π,π,π,
           60.0, 60.0,
           1]


res23_e1 = bbextract(bh_el1,
                     θ23,
                     lower23,
                     upper23;
                     TraceMode=:compact,
                     MaxFuncEvals=20000)
p = triptic(img_el1, res23_e1[1])
display(p)


res23_e2 = bbextract(bh_el2,
                     θ23,
                     lower23,
                     upper23;
                     TraceMode=:compact,
                     MaxFuncEvals=40000)
triptic(img_el2, res23_e2[1])

# Now to a 2/4 expansion
θ24 = CosineRing{2,4}(20.0,
                      [5.0,1.0],
                      [0.0],
                      0.2, 0.0,
                      [0.5, -0.5, 0.1,0.0],
                      [π/6, 0.0, 0.0,0.0],
                      -10.0, 10.0
                      ) + 1.0*Constant()

lower24 = [1.0,
           0.1, -5.0,
           -π,
           0.01, 0.0,
           0.01,-0.99,-0.99, -0.99,
           -π,-π,-π, -π,
           -60.0, -60.0,
           1e-6]

upper24 = [40.0,
           15.0, 5.0,
           π,
           0.99,π,
           0.99,0.99,0.99, 0.99,
           π,π,π,π,
           60.0, 60.0,
           1]


res24_e1 = bbextract(bh_el1,
                     θ24,
                     lower24,
                     upper24;
                     TraceMode=:compact,
                     MaxFuncEvals=20000)
triptic(img_el1, res24_e1[1])



res24_e2 = bbextract(bh_el2,
                     θ24,
                     lower24,
                     upper24;
                     TraceMode=:compact,
                     MaxFuncEvals=40000)
triptic(img_el2, res24_e2[1])



θ16 = CosineRing{1m,6}(20.0,
                  [5.0],
                  [],
                  0.2, 0.0,
                  [0.5, -0.5, 0.1, 0.1, 0.0, 0.0],
                  [π/6, 0.0, 0.0, 0.0, 0.0, 0.0],
                  -10.0, 10.0
                ) + 1.0*Constant()
lower16 = [1.0,
       0.1,

       0.01, 0.0,
       0.01,-0.99,-0.99,-0.99,-0.99,-0.99,
       -π,-π,-π,-π,-π,-π,
       -60.0, -60.0,
       1e-6]
upper16 = [40.0,
           10.0,
           0.99,π,
           0.99,0.99,0.99,0.99,0.99,0.99,
           π,π,π,π,π,π,
           60.0, 60.0,
           1]


res16_e1 = bbextract(bh_el1,
           θ16,
           lower16,
           upper16;
           TraceMode=:compact,
           MaxFuncEvals=20000)
triptic(img_el1, res16_e1[1])


res16_e2 = bbextract(bh_el2,
           θ16,
           lower16,
           upper16;
           TraceMode=:compact,
           MaxFuncEvals=20000)
triptic(img_el2, res16_e2[1])


#What about the GRMHD sims
resg_gr1 = bbextract(bh_gr1,
                    θg,
                    lowerg,
                    upperg;
                    TraceMode=:compact,
                    MaxFuncEvals=20000)
triptic(img_gr1, resg_gr1[1])

resg_gr2 = bbextract(bh_gr2,
                    θg,
                    lowerg,
                    upperg;
                    TraceMode=:compact,
                    MaxFuncEvals=20000)
triptic(img_gr2, resg_gr2[1])

#Well gr1 looks decent by 2 is a little off!
#gr 2 looks like a 2/4 expansion to my eye so lets try that


res24_gr2 = bbextract(bh_gr2,
                     θ24,
                     lower24,
                     upper24;
                     TraceMode=:compact,
                     MaxFuncEvals=40000)
triptic(img_gr2, res24_gr2[1])

#Lets also try the 14 expansion
res24_gr1 = bbextract(bh_gr1,
                     θ24,
                     lower24,
                     upper24;
                     TraceMode=:compact,
                     MaxFuncEvals=40000)
triptic(img_gr1, res24_gr1[1])

# MD5 sum of a file
md5sum(fn) = split(readstring(`md5sum $fn`))[1]

## tools.jl tests
#################
const todelete = []
function tempfn(ext="")
    out = tempname()*ext
    push!(todelete, out)
    out
end
function delfls()
    for f in todelete
        if isfile(f)
            rm(f)
        end
    end
end

## .agr readers:
g1 = BM._read_agr("testfiles/wiki.agr")
g1b = BM._read_agr("testfiles/wiki.bin")
@test g1.v == [-9999 -9999 5 2
                -9999 20 100 36
                3 8 35 10
                32 42 50 6
                88 75 27 9
                13 5 1 -9999]
@test g1.nc==4
@test g1.nr==6
@test g1.xll==0.0
@test g1.yll==-40.0
@test g1.dx==50.0
@test g1.NA==-9999
# write and read back:
tfn1 = tempfn(".agr")
write_agr(g1, tfn1)
g2 = BM._read_agr(tfn1, Float64)
@test g1==g2
tfn2 = tempfn(".agr")
write_agr(g2, tfn2)
@test md5sum(tfn1)==md5sum(tfn2)

# write and read back binary:
tfn3 = tempfn(".bin")
write_agr(g1, tfn3)
g2 = BM._read_agr(tfn3)
g1_ = BM._read_agr("testfiles/wiki.agr", Float32)
@test g1_==g2
@test md5sum(tfn3)==md5sum("testfiles/wiki.bin")

# larger file:
g3 = BM._read_agr("testfiles/t2.bin", Float32)
@test g3.nc==121
@test g3.nr==81
@test g3.xll==678750.0f0
@test g3.yll==140250.0f0
@test g3.dx==25.0f0
@test g3.NA==-9999.0f0

# larger file:
# @test_throws ErrorException BM._read_agr("testfiles/t2.bin", Float64)


tfn4 = tempfn(".bin")
write_agr(g3, tfn4)
@test md5sum("testfiles/t2.bin")==md5sum(tfn4)

# Gridded
gg = Gridded(g1)
@test BM.AGR(gg, NA_agr=-9999)==g1


## .xyn reader
po = read_xyn("testfiles/poly.xyn")
@test length(po)==1
@test size(po[1])==(2,24)

po = read_xyn("testfiles/multipoly.xyzn", hasz=true)
@test length(po)==4
@test size(po[1])==(3,4)
@test size(po[2])==(3,5)
@test size(po[3])==(3,4)
@test size(po[4])==(3,3)


# tidy up temp-files
delfls()


## IN poly
@test BM.leftorright(0.5,0.5, 1,0,1,1)==-1
@test BM.leftorright(1.5,.5, 1,0,1,1)==1
@test BM.leftorright(1,0.5, 1,0,1,1)==0

poc = BM.concat_poly(po)
@test !inpoly([0.0,0], poc[1])

poly = Float64[0 0
               0 1
               1 1
               1 0
               0 0]'
p1 = [0.5, 0.5]
p2 = [0.5, 0.99]
p22 = [0.5, 1] # on edge
p23 = [0.5, 0] # on edge
p24 = [0, 0]   # on corner
p25 = [0, .4]   # on edge
p3 = [0.5, 1.1]

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)

# clockwise poly
poly = Float64[0 0
               1 0
               1 1
               0 1
               0 0]'

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)


# cross-over poly
poly = Float64[0 0
               1 0
               0 1
               1 1
               0 0]'
if VERSION>=v"0.5-"
    eval(:(@test_broken inpoly(p1, poly) )) # should be true
end
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test !inpoly(p25, poly) # different
@test !inpoly(p3, poly)


# with interior region
poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside interior poly: i.e. labeled as outside
@test !inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.4 0.2
               0.2 0.2
               0.2 0.4
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.2 0.4
               0.2 0.2
               0.4 0.2
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior #1
               0.1 0.1
               0.1 0.6
               0.4 0.6
               0.4 0.6
               0.4 0.1
               0.1 0.1
               0 0
               # interior #2
               0.6 0.4
               0.6 0.6
               0.8 0.6
               0.8 0.4
               0.6 0.4
               0 0
               # exterior
               0 1
               1 1
               1 0
               0 0]'
@test !inpoly([0.2,0.4], poly)
@test !inpoly([0.3,0.15], poly)
@test inpoly([0.5,0.4], poly)
@test inpoly([0.5,0.2], poly)
@test !inpoly([0.7,0.5], poly)


#####
# boxcar
#####
#nr,nc = 4,5
T = Float32
nr,nc = 20,31
#nr,nc = 600,500
window = 3
# orig = (1.0:nr)''*(1:nc)'
# weights = ones(Int,nr,nc)
srand(1)
orig = rand(T,nr,nc)
weights = rand(nr,nc)
weightsb = bitrand(nr,nc)
weightsbb = convert(Matrix{Bool}, weightsb)
# with weights
filt1 = BM.boxcar(orig, window, weights)
M = BM.boxcar_matrix(T, window, weights)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = reshape(M*reshape(orig, length(orig)), size(orig))
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test_approx_eq filt1[i] filt2[i]
end
# with weightsb
filt1 = BM.boxcar(orig, window, weightsb)
M = BM.boxcar_matrix(T, window, weightsb)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = reshape(M*reshape(orig, length(orig)), size(orig))
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test_approx_eq filt1[i] filt2[i]
end
# with weightsbb
filt1 = BM.boxcar(orig, window, weightsbb)
M = BM.boxcar_matrix(T, window, weightsbb)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = reshape(M*reshape(orig, length(orig)), size(orig))
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test_approx_eq filt1[i] filt2[i]
end

####
# piecewiselinear
####
f = BM.piecewiselinear([78.0], [-12.0])
@test f(-1)==-12
@test f(100)==-12

f = BM.piecewiselinear([0, 1, 2], [-1, 1, -2])
@test f(0)==-1
@test f(0.5)==0
@test f(1)==1
@test f(1.5)==-0.5
@test f(2)==-2

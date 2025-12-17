using ITensors, ITensorMPS, CairoMakie, ColorSchemes, Printf, HDF5, DelimitedFiles

# Definitions of states and operators according to the {5}-representation

ITensors.space(::SiteType"site_basis") = 4;

ITensors.state(::StateName"0", ::SiteType"site_basis") = [1.0,0,0,0]
ITensors.state(::StateName"1", ::SiteType"site_basis") = [0,1.0,0,0]
ITensors.state(::StateName"2", ::SiteType"site_basis") = [0,0,1.0,0]
ITensors.state(::StateName"3", ::SiteType"site_basis") = [0,0,0,1.0]

ITensors.op(::OpName"0to1",::SiteType"site_basis")=
    [ 0  0   0   0
      1  0   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"0to2",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      1  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"0to3",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   0   0
      1  0   0   0 ]



ITensors.op(::OpName"1to0",::SiteType"site_basis")=
    [ 0  1   0   0
      0  0   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"1to2",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  1   0   0
      0  0   0   0 ]

ITensors.op(::OpName"1to3",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   0   0
      0  1   0   0 ]



ITensors.op(::OpName"2to0",::SiteType"site_basis")=
    [ 0  0   1   0
      0  0   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"2to1",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   1   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"2to3",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   0   0
      0  0   1   0 ]




ITensors.op(::OpName"3to0",::SiteType"site_basis")=
    [ 0  0   0   1
      0  0   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"3to1",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   1
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"3to2",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   0   1
      0  0   0   0 ]



ITensors.op(::OpName"p0",::SiteType"site_basis")=
    [ 1  0   0   0
      0  0   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"p1",::SiteType"site_basis")=
    [ 0  0   0   0
      0  1   0   0
      0  0   0   0
      0  0   0   0 ]

ITensors.op(::OpName"p2",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   1   0
      0  0   0   0 ]

ITensors.op(::OpName"p3",::SiteType"site_basis")=
    [ 0  0   0   0
      0  0   0   0
      0  0   0   0
      0  0   0   1 ]

ITensors.op(::OpName"num_is_flux",::SiteType"site_basis")=
    [ 0   0   0   0
      0   1   0   0
      0   0   1   0
      0   0   0   1 ]
  

ITensors.op(::OpName"right_f",::SiteType"site_basis")=
    [ 0   0   0   0
      0   1   0   0
      0   0   0   0
      0   0   0   1 ]

ITensors.op(::OpName"up_f",::SiteType"site_basis")=
    [ 0   0   0   0
      0   1   0   0
      0   0   1   0
      0   0   0   0 ]

ITensors.op(::OpName"down_f",::SiteType"site_basis")=
    [ 0   0   0   0
      0   0   0   0
      0   0   1   0
      0   0   0   1 ]


# definition of global variables

global Nx = parse(Int,ARGS[1]) # length in x-direction
global Ny = parse(Int,ARGS[2]) # length in y-direction from arguments
global N = Nx*Ny
global sites = siteinds("site_basis",N); # initialize site indices
global p = 10000 # constant for external charges penalty term

global T1 = -128 # constants {5}-representation hamiltonian
global T2 = 32
global T3 = 32
global T4 = 32
global T5 = -32
global T6 = -32
global T7 = -32
global T8 = 32

function get_H_mag() # generates magnetic hamiltonian MPO
    os = OpSum()
    for i=1:Ny-1
        for j=1:Nx-1 # applies the 64x64 plaquette hamiltonian for every plaquette
            os+=T1,"0to3",(i-1)*Nx+j,"0to2",i*Nx+j+1,"0to1",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"0to2",i*Nx+j+1,"1to0",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"0to2",i*Nx+j+1,"2to3",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"0to2",i*Nx+j+1,"3to2",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"1to3",i*Nx+j+1,"0to1",i*Nx+j
            os+=T6,"0to3",(i-1)*Nx+j,"1to3",i*Nx+j+1,"1to0",i*Nx+j
            os+=T7,"0to3",(i-1)*Nx+j,"1to3",i*Nx+j+1,"2to3",i*Nx+j
            os+=T5,"0to3",(i-1)*Nx+j,"1to3",i*Nx+j+1,"3to2",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"2to0",i*Nx+j+1,"0to1",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"2to0",i*Nx+j+1,"1to0",i*Nx+j
            os+=T6,"0to3",(i-1)*Nx+j,"2to0",i*Nx+j+1,"2to3",i*Nx+j
            os+=T4,"0to3",(i-1)*Nx+j,"2to0",i*Nx+j+1,"3to2",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"3to1",i*Nx+j+1,"0to1",i*Nx+j
            os+=T4,"0to3",(i-1)*Nx+j,"3to1",i*Nx+j+1,"1to0",i*Nx+j
            os+=T5,"0to3",(i-1)*Nx+j,"3to1",i*Nx+j+1,"2to3",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"3to1",i*Nx+j+1,"3to2",i*Nx+j
            os+=T2,"1to2",(i-1)*Nx+j,"0to2",i*Nx+j+1,"0to1",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"0to2",i*Nx+j+1,"1to0",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"0to2",i*Nx+j+1,"2to3",i*Nx+j
            os+=T7,"1to2",(i-1)*Nx+j,"0to2",i*Nx+j+1,"3to2",i*Nx+j
            os+=T3,"1to2",(i-1)*Nx+j,"1to3",i*Nx+j+1,"0to1",i*Nx+j
            os+=T7,"1to2",(i-1)*Nx+j,"1to3",i*Nx+j+1,"1to0",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"1to3",i*Nx+j+1,"2to3",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"1to3",i*Nx+j+1,"3to2",i*Nx+j
            os+=T4,"1to2",(i-1)*Nx+j,"2to0",i*Nx+j+1,"0to1",i*Nx+j
            os+=T2,"1to2",(i-1)*Nx+j,"2to0",i*Nx+j+1,"1to0",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"2to0",i*Nx+j+1,"2to3",i*Nx+j
            os+=T3,"1to2",(i-1)*Nx+j,"2to0",i*Nx+j+1,"3to2",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"3to1",i*Nx+j+1,"0to1",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"3to1",i*Nx+j+1,"1to0",i*Nx+j
            os+=T8,"1to2",(i-1)*Nx+j,"3to1",i*Nx+j+1,"2to3",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"3to1",i*Nx+j+1,"3to2",i*Nx+j
            os+=T2,"2to1",(i-1)*Nx+j,"0to2",i*Nx+j+1,"0to1",i*Nx+j
            os+=T4,"2to1",(i-1)*Nx+j,"0to2",i*Nx+j+1,"1to0",i*Nx+j
            os+=T3,"2to1",(i-1)*Nx+j,"0to2",i*Nx+j+1,"2to3",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"0to2",i*Nx+j+1,"3to2",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"1to3",i*Nx+j+1,"0to1",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"1to3",i*Nx+j+1,"1to0",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"1to3",i*Nx+j+1,"2to3",i*Nx+j
            os+=T8,"2to1",(i-1)*Nx+j,"1to3",i*Nx+j+1,"3to2",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"2to0",i*Nx+j+1,"0to1",i*Nx+j
            os+=T2,"2to1",(i-1)*Nx+j,"2to0",i*Nx+j+1,"1to0",i*Nx+j
            os+=T7,"2to1",(i-1)*Nx+j,"2to0",i*Nx+j+1,"2to3",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"2to0",i*Nx+j+1,"3to2",i*Nx+j
            os+=T7,"2to1",(i-1)*Nx+j,"3to1",i*Nx+j+1,"0to1",i*Nx+j
            os+=T3,"2to1",(i-1)*Nx+j,"3to1",i*Nx+j+1,"1to0",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"3to1",i*Nx+j+1,"2to3",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"3to1",i*Nx+j+1,"3to2",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"0to2",i*Nx+j+1,"0to1",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"0to2",i*Nx+j+1,"1to0",i*Nx+j
            os+=T4,"3to0",(i-1)*Nx+j,"0to2",i*Nx+j+1,"2to3",i*Nx+j
            os+=T6,"3to0",(i-1)*Nx+j,"0to2",i*Nx+j+1,"3to2",i*Nx+j
            os+=T4,"3to0",(i-1)*Nx+j,"1to3",i*Nx+j+1,"0to1",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"1to3",i*Nx+j+1,"1to0",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"1to3",i*Nx+j+1,"2to3",i*Nx+j
            os+=T5,"3to0",(i-1)*Nx+j,"1to3",i*Nx+j+1,"3to2",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"2to0",i*Nx+j+1,"0to1",i*Nx+j
            os+=T1,"3to0",(i-1)*Nx+j,"2to0",i*Nx+j+1,"1to0",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"2to0",i*Nx+j+1,"2to3",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"2to0",i*Nx+j+1,"3to2",i*Nx+j
            os+=T6,"3to0",(i-1)*Nx+j,"3to1",i*Nx+j+1,"0to1",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"3to1",i*Nx+j+1,"1to0",i*Nx+j
            os+=T5,"3to0",(i-1)*Nx+j,"3to1",i*Nx+j+1,"2to3",i*Nx+j
            os+=T7,"3to0",(i-1)*Nx+j,"3to1",i*Nx+j+1,"3to2",i*Nx+j
        end
    end
    for i=1:Ny-1
        j=Nx # additional plaquettes due to periodic boundries
            os+=T1,"0to3",(i-1)*Nx+j,"0to2",i*Nx+1,"0to1",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"0to2",i*Nx+1,"1to0",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"0to2",i*Nx+1,"2to3",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"0to2",i*Nx+1,"3to2",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"1to3",i*Nx+1,"0to1",i*Nx+j
            os+=T6,"0to3",(i-1)*Nx+j,"1to3",i*Nx+1,"1to0",i*Nx+j
            os+=T7,"0to3",(i-1)*Nx+j,"1to3",i*Nx+1,"2to3",i*Nx+j
            os+=T5,"0to3",(i-1)*Nx+j,"1to3",i*Nx+1,"3to2",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"2to0",i*Nx+1,"0to1",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"2to0",i*Nx+1,"1to0",i*Nx+j
            os+=T6,"0to3",(i-1)*Nx+j,"2to0",i*Nx+1,"2to3",i*Nx+j
            os+=T4,"0to3",(i-1)*Nx+j,"2to0",i*Nx+1,"3to2",i*Nx+j
            os+=T2,"0to3",(i-1)*Nx+j,"3to1",i*Nx+1,"0to1",i*Nx+j
            os+=T4,"0to3",(i-1)*Nx+j,"3to1",i*Nx+1,"1to0",i*Nx+j
            os+=T5,"0to3",(i-1)*Nx+j,"3to1",i*Nx+1,"2to3",i*Nx+j
            os+=T3,"0to3",(i-1)*Nx+j,"3to1",i*Nx+1,"3to2",i*Nx+j
            os+=T2,"1to2",(i-1)*Nx+j,"0to2",i*Nx+1,"0to1",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"0to2",i*Nx+1,"1to0",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"0to2",i*Nx+1,"2to3",i*Nx+j
            os+=T7,"1to2",(i-1)*Nx+j,"0to2",i*Nx+1,"3to2",i*Nx+j
            os+=T3,"1to2",(i-1)*Nx+j,"1to3",i*Nx+1,"0to1",i*Nx+j
            os+=T7,"1to2",(i-1)*Nx+j,"1to3",i*Nx+1,"1to0",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"1to3",i*Nx+1,"2to3",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"1to3",i*Nx+1,"3to2",i*Nx+j
            os+=T4,"1to2",(i-1)*Nx+j,"2to0",i*Nx+1,"0to1",i*Nx+j
            os+=T2,"1to2",(i-1)*Nx+j,"2to0",i*Nx+1,"1to0",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"2to0",i*Nx+1,"2to3",i*Nx+j
            os+=T3,"1to2",(i-1)*Nx+j,"2to0",i*Nx+1,"3to2",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"3to1",i*Nx+1,"0to1",i*Nx+j
            os+=T5,"1to2",(i-1)*Nx+j,"3to1",i*Nx+1,"1to0",i*Nx+j
            os+=T8,"1to2",(i-1)*Nx+j,"3to1",i*Nx+1,"2to3",i*Nx+j
            os+=T6,"1to2",(i-1)*Nx+j,"3to1",i*Nx+1,"3to2",i*Nx+j
            os+=T2,"2to1",(i-1)*Nx+j,"0to2",i*Nx+1,"0to1",i*Nx+j
            os+=T4,"2to1",(i-1)*Nx+j,"0to2",i*Nx+1,"1to0",i*Nx+j
            os+=T3,"2to1",(i-1)*Nx+j,"0to2",i*Nx+1,"2to3",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"0to2",i*Nx+1,"3to2",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"1to3",i*Nx+1,"0to1",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"1to3",i*Nx+1,"1to0",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"1to3",i*Nx+1,"2to3",i*Nx+j
            os+=T8,"2to1",(i-1)*Nx+j,"1to3",i*Nx+1,"3to2",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"2to0",i*Nx+1,"0to1",i*Nx+j
            os+=T2,"2to1",(i-1)*Nx+j,"2to0",i*Nx+1,"1to0",i*Nx+j
            os+=T7,"2to1",(i-1)*Nx+j,"2to0",i*Nx+1,"2to3",i*Nx+j
            os+=T5,"2to1",(i-1)*Nx+j,"2to0",i*Nx+1,"3to2",i*Nx+j
            os+=T7,"2to1",(i-1)*Nx+j,"3to1",i*Nx+1,"0to1",i*Nx+j
            os+=T3,"2to1",(i-1)*Nx+j,"3to1",i*Nx+1,"1to0",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"3to1",i*Nx+1,"2to3",i*Nx+j
            os+=T6,"2to1",(i-1)*Nx+j,"3to1",i*Nx+1,"3to2",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"0to2",i*Nx+1,"0to1",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"0to2",i*Nx+1,"1to0",i*Nx+j
            os+=T4,"3to0",(i-1)*Nx+j,"0to2",i*Nx+1,"2to3",i*Nx+j
            os+=T6,"3to0",(i-1)*Nx+j,"0to2",i*Nx+1,"3to2",i*Nx+j
            os+=T4,"3to0",(i-1)*Nx+j,"1to3",i*Nx+1,"0to1",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"1to3",i*Nx+1,"1to0",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"1to3",i*Nx+1,"2to3",i*Nx+j
            os+=T5,"3to0",(i-1)*Nx+j,"1to3",i*Nx+1,"3to2",i*Nx+j
            os+=T3,"3to0",(i-1)*Nx+j,"2to0",i*Nx+1,"0to1",i*Nx+j
            os+=T1,"3to0",(i-1)*Nx+j,"2to0",i*Nx+1,"1to0",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"2to0",i*Nx+1,"2to3",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"2to0",i*Nx+1,"3to2",i*Nx+j
            os+=T6,"3to0",(i-1)*Nx+j,"3to1",i*Nx+1,"0to1",i*Nx+j
            os+=T2,"3to0",(i-1)*Nx+j,"3to1",i*Nx+1,"1to0",i*Nx+j
            os+=T5,"3to0",(i-1)*Nx+j,"3to1",i*Nx+1,"2to3",i*Nx+j
            os+=T7,"3to0",(i-1)*Nx+j,"3to1",i*Nx+1,"3to2",i*Nx+j
    end
    H_mag = MPO(os,sites)
    return H_mag
end

function get_H_el() # generates electric hamiltonian MPO
    os_el = OpSum()
    for i=1:N
        os_el+=6/4*2,"num_is_flux",i # the square of L_i and R_i have 1/4 eig.val if the link carries flux, else 0. There are 2 flux carrying links in states 1, 2 and 3
    end
    H_el = MPO(os_el,sites)
    return H_el
end

function get_H_gauss() # generates gauss's law penatly term
    os_gauss = OpSum()
    for i=1:Ny-1
        for j=1:Nx-1 # assigns EV of 1 to each prohibited configuration, EV of 0 else.
            os_gauss+=0,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+j+1,"p3",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p0",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p1",i*Nx+j+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p2",i*Nx+j+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+j+1,"p3",i*Nx+j+1
        end
    end
    for i=1:Ny-1
        j=Nx # additional Gauss's conditions from periodic boundries
            os_gauss+=0,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p0",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p1",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p2",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p0",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p1",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p2",(i-1)*Nx+1,"p3",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p0",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p1",i*Nx+1
            os_gauss+=1,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p2",i*Nx+1
            os_gauss+=0,"p3",(i-1)*Nx+j,"p3",(i-1)*Nx+1,"p3",i*Nx+1
    end
    H_gauss = MPO(os_gauss,sites)
    return H_gauss
end

function get_H_gauss_closed1() # generates modified Gauss's law for closed boundry at the bottom
    os_gauss_closed1 = OpSum()
    i=Ny-1
    for j=1:Nx-1 # assigns EV of 1 to each prohibited configuration, EV of 0 else.
        os_gauss_closed1+=0,"p0",i*Nx+j,"p0",i*Nx+j+1
        os_gauss_closed1+=0,"p0",i*Nx+j,"p1",i*Nx+j+1
        os_gauss_closed1+=1,"p0",i*Nx+j,"p2",i*Nx+j+1
        os_gauss_closed1+=1,"p0",i*Nx+j,"p3",i*Nx+j+1
    
        os_gauss_closed1+=1,"p1",i*Nx+j,"p0",i*Nx+j+1
        os_gauss_closed1+=1,"p1",i*Nx+j,"p1",i*Nx+j+1
        os_gauss_closed1+=0,"p1",i*Nx+j,"p2",i*Nx+j+1
        os_gauss_closed1+=0,"p1",i*Nx+j,"p3",i*Nx+j+1
    
        os_gauss_closed1+=0,"p2",i*Nx+j,"p0",i*Nx+j+1
        os_gauss_closed1+=0,"p2",i*Nx+j,"p1",i*Nx+j+1
        os_gauss_closed1+=1,"p2",i*Nx+j,"p2",i*Nx+j+1
        os_gauss_closed1+=1,"p2",i*Nx+j,"p3",i*Nx+j+1
    
        os_gauss_closed1+=1,"p3",i*Nx+j,"p0",i*Nx+j+1
        os_gauss_closed1+=1,"p3",i*Nx+j,"p1",i*Nx+j+1
        os_gauss_closed1+=0,"p3",i*Nx+j,"p2",i*Nx+j+1
        os_gauss_closed1+=0,"p3",i*Nx+j,"p3",i*Nx+j+1
    end

    i=Ny-1
    j=Nx # additional conditions from periodic boundries
        os_gauss_closed1+=0,"p0",i*Nx+j,"p0",i*Nx+1
        os_gauss_closed1+=0,"p0",i*Nx+j,"p1",i*Nx+1
        os_gauss_closed1+=1,"p0",i*Nx+j,"p2",i*Nx+1
        os_gauss_closed1+=1,"p0",i*Nx+j,"p3",i*Nx+1
    
        os_gauss_closed1+=1,"p1",i*Nx+j,"p0",i*Nx+1
        os_gauss_closed1+=1,"p1",i*Nx+j,"p1",i*Nx+1
        os_gauss_closed1+=0,"p1",i*Nx+j,"p2",i*Nx+1
        os_gauss_closed1+=0,"p1",i*Nx+j,"p3",i*Nx+1
    
        os_gauss_closed1+=0,"p2",i*Nx+j,"p0",i*Nx+1
        os_gauss_closed1+=0,"p2",i*Nx+j,"p1",i*Nx+1
        os_gauss_closed1+=1,"p2",i*Nx+j,"p2",i*Nx+1
        os_gauss_closed1+=1,"p2",i*Nx+j,"p3",i*Nx+1
    
        os_gauss_closed1+=1,"p3",i*Nx+j,"p0",i*Nx+1
        os_gauss_closed1+=1,"p3",i*Nx+j,"p1",i*Nx+1
        os_gauss_closed1+=0,"p3",i*Nx+j,"p2",i*Nx+1
        os_gauss_closed1+=0,"p3",i*Nx+j,"p3",i*Nx+1

    H_gauss_closed1 = MPO(os_gauss_closed1,sites)
    return H_gauss_closed1
end

function get_H_gauss_closed2() # generates penalty term for closed boundry at the top
    os_gauss_closed2 = OpSum()

    for j=1:Nx
        os_gauss_closed2+=0,"p0",j
        os_gauss_closed2+=1,"p1",j
        os_gauss_closed2+=1,"p2",j
        os_gauss_closed2+=0,"p3",j
    end

    H_gauss_closed2 = MPO(os_gauss_closed2,sites)
    return H_gauss_closed2
end

function charge_half_bottom(m) #add to closed1 to attach 1/2 charge on m-th link at the top
    i=Ny-1
    j = m
    os_gauss_closed3 = OpSum()
    os_gauss_closed3+=p,"p0",i*Nx+j,"p0",i*Nx+j+1
    os_gauss_closed3+=p,"p0",i*Nx+j,"p1",i*Nx+j+1
    os_gauss_closed3+=-1,"p0",i*Nx+j,"p2",i*Nx+j+1
    os_gauss_closed3+=-1,"p0",i*Nx+j,"p3",i*Nx+j+1
    os_gauss_closed3+=-1,"p1",i*Nx+j,"p0",i*Nx+j+1
    os_gauss_closed3+=-1,"p1",i*Nx+j,"p1",i*Nx+j+1
    os_gauss_closed3+=p,"p1",i*Nx+j,"p2",i*Nx+j+1
    os_gauss_closed3+=p,"p1",i*Nx+j,"p3",i*Nx+j+1
    os_gauss_closed3+=p,"p2",i*Nx+j,"p0",i*Nx+j+1
    os_gauss_closed3+=p,"p2",i*Nx+j,"p1",i*Nx+j+1
    os_gauss_closed3+=-1,"p2",i*Nx+j,"p2",i*Nx+j+1
    os_gauss_closed3+=-1,"p2",i*Nx+j,"p3",i*Nx+j+1
    os_gauss_closed3+=-1,"p3",i*Nx+j,"p0",i*Nx+j+1
    os_gauss_closed3+=-1,"p3",i*Nx+j,"p1",i*Nx+j+1
    os_gauss_closed3+=p,"p3",i*Nx+j,"p2",i*Nx+j+1
    os_gauss_closed3+=p,"p3",i*Nx+j,"p3",i*Nx+j+1
    return MPO(os_gauss_closed3,sites)
end

function charge_half_top(m) #add to closed2 to attach 1/2 charge on m-th link at the bottom
    j=m
    os_gauss_closed4 = OpSum()
    os_gauss_closed4+=p,"p0",j
    os_gauss_closed4+=-1,"p1",j
    os_gauss_closed4+=-1,"p2",j
    os_gauss_closed4+=p,"p3",j
    return MPO(os_gauss_closed4,sites)
end

function modup(x) # helping tool to compensate for the fact that julia starts arrays on 1 and thus indexing with mod() is a bit of a pain
    return mod(x-1,Nx)+1
end

function get_groundstate(g2,K,H_el,H_mag,H_gauss,H_gauss_closed1,H_gauss_closed2, charge_bot, charge_top,linkDim=max(20, ceil(Int,Nx*Ny*1.25)); rand_init=false, minsweeps=50, noise_sweeps=15, pass_init=false, init_state=0) # performs DMRG for given g². Returns energy and state
    if rand_init # option to start with a random initial state
        init = random_mps(sites;linkdims=linkDim)
    elseif pass_init # option to pass an initial state
        init = init_state
        linkDim = max(maxlinkdim(init),linkDim) # adjust initial link dimension based on link dimension of passed state
        linkDim = min(linkDim, 256) # however at most 256
    else # generates initial state of 50% vaccuum state and 50% string state between the two external charges
        if Ny%2 == 0
            init1 = MPS(sites,["0" for n in 1:N])

            vec_init = ["0" for n in 1:N]
            for i in 0:Ny-1
                if i%2 == 0
                    vec_init[i*Nx+modup(3-trunc(Int,Ny/2)+fld(i,2))] = "2"
                else
                    vec_init[i*Nx+modup(3-trunc(Int,Ny/2)+fld(i,2))] = "1"
                end
            end

            init2 = MPS(sites,vec_init)

            init = 0.5*init1+0.5*init2
        else
            init1 = MPS(sites,["0" for n in 1:N])

            vec_init = ["0" for n in 1:N]
            vec_init[modup(3-trunc(Int,Ny/2))] = "2"
            for i in 0:Ny-2
                if i%2 == 0
                    vec_init[(i+1)*Nx+modup(3-trunc(Int,Ny/2)+fld(i,2))] = "2"
                else
                    vec_init[(i+1)*Nx+modup(3-trunc(Int,Ny/2)+fld(i,2))] = "1"
                end
            end

            init2 = MPS(sites,vec_init)

            init = 0.5*init1+0.5*init2
        end
    end
    
    obs = DMRGObserver(energy_tol=5E-5, minsweeps=minsweeps) # terminate DMRG when energy changes less than 5E-5 and at least minsweep (default is 50) sweeps have been performed
    maxdim = 10000 # limit maximum linkdim to avoid running out of ram
    energy,psi = dmrg([g2*H_el, 1/g2*H_mag,K*H_gauss,K*(H_gauss_closed1+charge_bot),K*(H_gauss_closed2+charge_top)], init ;nsweeps=10000,cutoff=1E-9, mindim=5, eigsolve_tol=1e-14,maxdim=Int.(reduce(vcat,[
    fill(min(maxdim,8*linkDim/8),3), # define progressively growing maximum link dimension
    fill(min(maxdim,8*linkDim/4),3),
    fill(min(maxdim,8*linkDim),3),
    fill(min(maxdim,8*linkDim*2),3),
    fill(min(maxdim,8*linkDim*3),3),
    fill(min(maxdim,8*linkDim*5),3),
    fill(min(maxdim,8*linkDim*10),5),
    fill(min(maxdim,8*linkDim*20),5),
    fill(min(maxdim,8*linkDim*50),5),
    fill(min(maxdim,8*linkDim*100),10),
    fill(min(maxdim,8*linkDim*200),10),
    fill(min(maxdim,8*linkDim*500),10),
    ],)), noise=reduce(vcat, [fill(1E-4,noise_sweeps), [0.0]]),observer=obs)

    return psi, energy
end

function color(i) # colors for plotting
    cmap2 = to_colormap(Reverse(:thermal))
    return Makie.interpolated_getindex(cmap2, i, (0,1))
end


function save_graph_links(state, name, title) # generates heatmap of expectation values of projector 1 - |0> <0|
    sites = []
    for i=1:Ny
        for j=1:Nx
            push!(sites,(sqrt(3)*(j-(i-1)/2-1),-1.5*(i-1)))
        end
    end
    p_sites = Point2f[sites...]

    right = (sqrt(3)/2,-0.5)
    up = (0.0,1.0)
    down = (-sqrt(3)/2,-0.5)

    right_line = fill(right,N)
    up_line = fill(up,N)
    down_line = fill(down,N)
    for i=1:N
        right_line[i] = right_line[i] .+ sites[i]
        up_line[i] = up_line[i] .+ sites[i]
        down_line[i] = down_line[i] .+ sites[i]
    end

    p_right = Point2f[right_line...]
    p_up = Point2f[up_line...]
    p_down = Point2f[down_line...]

    f_right = expect(state,"right_f")
    c_right = [color(i) for i in f_right]

    f_up = expect(state,"up_f")
    c_up = [color(i) for i in f_up]

    f_down = expect(state,"down_f")
    c_down = [color(i) for i in f_down]

    f = Figure(size = (1.5*80*(Nx+Ny/2), 1.5*69*Ny))

    ax = Axis(f[1,1], autolimitaspect=1)

    hidedecorations!(ax)


    for i=1:N
        lines!(stack([p_sites[i],p_right[i]]), color=c_right[i],linewidth=10)
        lines!(stack([p_sites[i],p_up[i]]), color=c_up[i],linewidth=10)
        lines!(stack([p_sites[i],p_down[i]]), color=c_down[i],linewidth=10)
    end

    Colorbar(f[1,2], limits=(0,1), colormap=Reverse(:thermal))

    scatter!(p_sites, color="black", markersize=19, marker=:dtriangle)

    save(name*".pdf", f, pdf_version="1.4")

end

function ground_sweep(g2s; carryover=false, w=1)
    l = length(g2s)
    # generate all hamiltonians and penalty terms
    H_el = get_H_el()
    H_mag = get_H_mag()
    H_gauss_closed1 = get_H_gauss_closed1()
    H_gauss_closed2 = get_H_gauss_closed2()
    H_gauss = get_H_gauss()
    charge_bot = charge_half_bottom(2)
    charge_top = charge_half_top(modup(3-floor(Int,Ny/2))) # adjust charge position such that string runs perpendicular to y-boundry

    Ks = 100 .*[max(g2,1/g2) for g2 in g2s] # set constant for penalty terms based on parameter regime

    # initialize all the output lists
    gauss=[]
    closed1=[]
    closed2=[]
    variance=[]

    ex0s = []
    ex1s = []
    ex2s = []
    ex3s = []

    Es = []
    init = 0

    println("Using a "*string(Nx)*" by "*string(Ny)*" lattice:")

    for i=w:l
        println(string(i)*"/"*string(l))
        if carryover
            if i == w # first iteration has minsweeps of 50 and starts with prepared string state
                state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed1,H_gauss_closed2, charge_bot, charge_top ,noise_sweeps=15,rand_init=false)
                init = state
            else # any further iteration has minsweeps of 7 and takes previous state as initial state
                state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed1,H_gauss_closed2, charge_bot, charge_top, pass_init=true, init_state=init ,noise_sweeps=4, minsweeps=7)
                init = state
            end
        else # option for just independently using DMRG with random initial state
            state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed1,H_gauss_closed2, charge_bot, charge_top ,noise_sweeps=25)
        end

        # calculate and save observables

        push!(gauss,inner(state',H_gauss,state))
        push!(closed1,inner(state',H_gauss_closed1+charge_bot,state))
        push!(closed2,inner(state',H_gauss_closed2+charge_top,state))
        push!(variance,inner(g2s[i]*H_el+1/g2s[i]*H_mag,state,g2s[i]*H_el+1/g2s[i]*H_mag,state)-energy^2)
        push!(Es,energy)
        push!(ex0s,expect(state,"p0"))
        push!(ex1s,expect(state,"p1"))
        push!(ex2s,expect(state,"p2"))
        push!(ex3s,expect(state,"p3"))

        save_graph_links(state,"link_"*string(Nx)*"x"*string(Ny)*"_g2_"*@sprintf("%2.2f", g2s[i]), @sprintf("g=%.2f", g2s[i]))

        writedlm("AAA_g_"*string(Nx)*"x"*string(Ny)*".txt",g2s)
        writedlm("AAA_E_"*string(Nx)*"x"*string(Ny)*".txt",Es)
        writedlm("AAA_gauss_"*string(Nx)*"x"*string(Ny)*".txt",gauss)
        writedlm("AAA_closed1_"*string(Nx)*"x"*string(Ny)*".txt",closed1)
        writedlm("AAA_closed2_"*string(Nx)*"x"*string(Ny)*".txt",closed2)
        writedlm("AAA_variance_"*string(Nx)*"x"*string(Ny)*".txt",variance)
        writedlm("AAA_ex0_"*string(Nx)*"x"*string(Ny)*".txt", ex0s)
        writedlm("AAA_ex1_"*string(Nx)*"x"*string(Ny)*".txt", ex1s)
        writedlm("AAA_ex2_"*string(Nx)*"x"*string(Ny)*".txt", ex2s)
        writedlm("AAA_ex3_"*string(Nx)*"x"*string(Ny)*".txt", ex3s)

    end

end


try
    ground_sweep(LinRange(8,0.5,31),carryover=true, w=parse(Int,ARGS[3])) # call julia dmrg_5xN.jl N_y w to start at the w-th g2 value. Usefull after crashes or timeouts
catch e
    ground_sweep(LinRange(8,0.5,31),carryover=true) # normal call that just runs through in order
end

#ground_sweep(LinRange(parse(Float64,ARGS[3]),parse(Float64,ARGS[3]),1),carryover=true) # call julia dmrg_5xN.jl
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
      0   0   0   0
      0   0   1   0
      0   0   0   1 ]

ITensors.op(::OpName"up_f",::SiteType"site_basis")=
    [ 0   0   0   0
      0   1   0   0
      0   0   1   0
      0   0   0   0 ]

ITensors.op(::OpName"down_f",::SiteType"site_basis")=
    [ 0   0   0   0
      0   1   0   0
      0   0   0   0
      0   0   0   1 ]


# definition of global variables

global Nx = parse(Int,ARGS[1]) # length in x-direction (must be even!)
global Ny = parse(Int,ARGS[2]) # length in y-direction from arguments
global N = Nx*Ny
@assert Nx%2==0 "Nx must be even"
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

function get_OS_plaquette(os, x, y)
    i = y-1
    if x%2==0 && x!=Nx # x!=Nx takes care of periodic boundary
        i2 = i+1
    else
        i2 = i
    end
    j = x

    os+=T1,"0to3",i*Nx+j,"0to1",i2*Nx+j+1,"0to2",(i+1)*Nx+j

    os+=T1,"3to0",i*Nx+j,"1to0",i2*Nx+j+1,"2to0",(i+1)*Nx+j

    os+=T2,"2to1",i*Nx+j,"0to1",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T2,"0to3",i*Nx+j,"3to2",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T2,"0to3",i*Nx+j,"0to1",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T2,"1to2",i*Nx+j,"0to1",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T2,"0to3",i*Nx+j,"2to3",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T2,"0to3",i*Nx+j,"0to1",i2*Nx+j+1,"3to1",(i+1)*Nx+j

    os+=T2,"1to2",i*Nx+j,"1to0",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T2,"3to0",i*Nx+j,"2to3",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T2,"3to0",i*Nx+j,"1to0",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T2,"2to1",i*Nx+j,"1to0",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T2,"3to0",i*Nx+j,"3to2",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T2,"3to0",i*Nx+j,"1to0",i2*Nx+j+1,"1to3",(i+1)*Nx+j

    os+=T3,"2to1",i*Nx+j,"2to3",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T3,"0to3",i*Nx+j,"3to2",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T3,"1to2",i*Nx+j,"0to1",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T3,"3to0",i*Nx+j,"0to1",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T3,"0to3",i*Nx+j,"1to0",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T3,"0to3",i*Nx+j,"0to1",i2*Nx+j+1,"2to0",(i+1)*Nx+j

    os+=T3,"1to2",i*Nx+j,"3to2",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T3,"3to0",i*Nx+j,"2to3",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T3,"2to1",i*Nx+j,"1to0",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T3,"0to3",i*Nx+j,"1to0",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T3,"3to0",i*Nx+j,"0to1",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T3,"3to0",i*Nx+j,"1to0",i2*Nx+j+1,"0to2",(i+1)*Nx+j

    os+=T4,"2to1",i*Nx+j,"1to0",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T4,"0to3",i*Nx+j,"3to2",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T4,"3to0",i*Nx+j,"0to1",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    
    os+=T4,"1to2",i*Nx+j,"0to1",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T4,"3to0",i*Nx+j,"2to3",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T4,"0to3",i*Nx+j,"1to0",i2*Nx+j+1,"3to1",(i+1)*Nx+j

    os+=T5,"2to1",i*Nx+j,"3to2",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T5,"0to3",i*Nx+j,"3to2",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T5,"2to1",i*Nx+j,"0to1",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T5,"1to2",i*Nx+j,"2to3",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T5,"0to3",i*Nx+j,"2to3",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T5,"1to2",i*Nx+j,"0to1",i2*Nx+j+1,"3to1",(i+1)*Nx+j

    os+=T5,"1to2",i*Nx+j,"2to3",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T5,"3to0",i*Nx+j,"2to3",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T5,"1to2",i*Nx+j,"1to0",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T5,"2to1",i*Nx+j,"3to2",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T5,"3to0",i*Nx+j,"3to2",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T5,"2to1",i*Nx+j,"1to0",i2*Nx+j+1,"1to3",(i+1)*Nx+j

    os+=T6,"2to1",i*Nx+j,"3to2",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T6,"1to2",i*Nx+j,"3to2",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T6,"2to1",i*Nx+j,"2to3",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T6,"0to3",i*Nx+j,"2to3",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T6,"3to0",i*Nx+j,"0to1",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T6,"1to2",i*Nx+j,"1to0",i2*Nx+j+1,"0to2",(i+1)*Nx+j

    os+=T6,"1to2",i*Nx+j,"2to3",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T6,"2to1",i*Nx+j,"2to3",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T6,"1to2",i*Nx+j,"3to2",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T6,"3to0",i*Nx+j,"3to2",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T6,"0to3",i*Nx+j,"1to0",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T6,"2to1",i*Nx+j,"0to1",i2*Nx+j+1,"2to0",(i+1)*Nx+j

    os+=T7,"2to1",i*Nx+j,"0to1",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    os+=T7,"1to2",i*Nx+j,"3to2",i2*Nx+j+1,"0to2",(i+1)*Nx+j
    os+=T7,"0to3",i*Nx+j,"2to3",i2*Nx+j+1,"1to3",(i+1)*Nx+j

    os+=T7,"1to2",i*Nx+j,"1to0",i2*Nx+j+1,"1to3",(i+1)*Nx+j
    os+=T7,"2to1",i*Nx+j,"2to3",i2*Nx+j+1,"2to0",(i+1)*Nx+j
    os+=T7,"3to0",i*Nx+j,"3to2",i2*Nx+j+1,"3to1",(i+1)*Nx+j

    os+=T8,"2to1",i*Nx+j,"3to2",i2*Nx+j+1,"1to3",(i+1)*Nx+j

    os+=T8,"1to2",i*Nx+j,"2to3",i2*Nx+j+1,"3to1",(i+1)*Nx+j
    return os
end

function get_H_mag() # generates magnetic hamiltonian MPO
    os = OpSum()
    for y=1:Ny-1
        for x=1:Nx # applies the 64x64 plaquette hamiltonian to every plaquette
            os=get_OS_plaquette(os, x, y)
        end
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

function get_OS_gauss_site(os, x, y)
    i = y-1
    if x%2==0 && x!=Nx # x!=Nx takes care of periodic boundary
        i2 = i
        i3 = i+1
    else
        i2 = i-1
        i3 = i
    end
    j = x

    os+=0,"p0",i*Nx+j,"p0",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p0",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p0",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p0",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p1",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p1",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p1",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p1",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p2",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p2",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p2",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p2",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p3",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p3",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p0",i*Nx+j,"p3",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p0",i*Nx+j,"p3",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p0",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p0",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p0",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p0",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p1",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p1",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p1",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p1",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p2",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p2",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p2",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p2",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p3",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p3",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p1",i*Nx+j,"p3",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p1",i*Nx+j,"p3",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p0",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p0",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p0",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p0",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p1",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p1",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p1",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p1",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p2",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p2",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p2",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p2",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p3",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p3",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p2",i*Nx+j,"p3",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p2",i*Nx+j,"p3",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p0",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p0",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p0",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p0",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p1",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p1",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p1",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p1",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p2",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p2",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p2",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p2",i2*Nx+j+1,"p3",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p3",i2*Nx+j+1,"p0",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p3",i2*Nx+j+1,"p1",i3*Nx+j+1
    os+=1,"p3",i*Nx+j,"p3",i2*Nx+j+1,"p2",i3*Nx+j+1
    os+=0,"p3",i*Nx+j,"p3",i2*Nx+j+1,"p3",i3*Nx+j+1
    return os
end

function get_H_gauss() # generates gauss's law penatly term
    os = OpSum()
    for y=2:Ny-1
        for x=1:Nx # assigns EV of 1 to each prohibited configuration, EV of 0 else.
            os=get_OS_gauss_site(os, x, y)
        end
    end
    for x=2:2:Nx
        y=1
        os=get_OS_gauss_site(os, x, y)
    end
    for x=1:2:Nx
        y=Ny
        os=get_OS_gauss_site(os, x, y)
    end
    H_gauss = MPO(os,sites)
    return H_gauss
end

function get_H_gauss_closed_bottom() # generates modified Gauss's law for closed boundary at the bottom
    # assigns EV of 1 to each prohibited configuration, EV of 0 else.
    os = OpSum()
    i=Ny-1
    for x=2:2:Nx
        j=x
        if x==Nx
            j2=1
        else
            j2=x+1
        end
        os+=0,"p0",i*Nx+j
        os+=1,"p1",i*Nx+j
        os+=0,"p2",i*Nx+j
        os+=1,"p3",i*Nx+j

        os+=0,"p0",i*Nx+j,"p0",i*Nx+j2
        os+=1,"p0",i*Nx+j,"p1",i*Nx+j2
        os+=0,"p0",i*Nx+j,"p2",i*Nx+j2
        os+=1,"p0",i*Nx+j,"p3",i*Nx+j2
        os+=0,"p1",i*Nx+j,"p0",i*Nx+j2
        os+=1,"p1",i*Nx+j,"p1",i*Nx+j2
        os+=0,"p1",i*Nx+j,"p2",i*Nx+j2
        os+=1,"p1",i*Nx+j,"p3",i*Nx+j2
        os+=1,"p2",i*Nx+j,"p0",i*Nx+j2
        os+=0,"p2",i*Nx+j,"p1",i*Nx+j2
        os+=1,"p2",i*Nx+j,"p2",i*Nx+j2
        os+=0,"p2",i*Nx+j,"p3",i*Nx+j2
        os+=1,"p3",i*Nx+j,"p0",i*Nx+j2
        os+=0,"p3",i*Nx+j,"p1",i*Nx+j2
        os+=1,"p3",i*Nx+j,"p2",i*Nx+j2
        os+=0,"p3",i*Nx+j,"p3",i*Nx+j2
    end
    H_gauss_closed_bottom = MPO(os,sites)
    return H_gauss_closed_bottom
end

function get_H_gauss_closed_top() # generates modified Gauss's law for closed boundary at the top
    # assigns EV of 1 to each prohibited configuration, EV of 0 else.
    os = OpSum()
    i=0
    for x=1:2:Nx
        j=x
        os+=0,"p0",i*Nx+j
        os+=1,"p1",i*Nx+j
        os+=1,"p2",i*Nx+j
        os+=0,"p3",i*Nx+j

        os+=0,"p0",i*Nx+j,"p0",i*Nx+j+1
        os+=1,"p0",i*Nx+j,"p1",i*Nx+j+1
        os+=1,"p0",i*Nx+j,"p2",i*Nx+j+1
        os+=0,"p0",i*Nx+j,"p3",i*Nx+j+1
        os+=0,"p1",i*Nx+j,"p0",i*Nx+j+1
        os+=1,"p1",i*Nx+j,"p1",i*Nx+j+1
        os+=1,"p1",i*Nx+j,"p2",i*Nx+j+1
        os+=0,"p1",i*Nx+j,"p3",i*Nx+j+1
        os+=1,"p2",i*Nx+j,"p0",i*Nx+j+1
        os+=0,"p2",i*Nx+j,"p1",i*Nx+j+1
        os+=0,"p2",i*Nx+j,"p2",i*Nx+j+1
        os+=1,"p2",i*Nx+j,"p3",i*Nx+j+1
        os+=1,"p3",i*Nx+j,"p0",i*Nx+j+1
        os+=0,"p3",i*Nx+j,"p1",i*Nx+j+1
        os+=0,"p3",i*Nx+j,"p2",i*Nx+j+1
        os+=1,"p3",i*Nx+j,"p3",i*Nx+j+1
    end
    H_gauss_closed_top = MPO(os,sites)
    return H_gauss_closed_top
end

function get_groundstate(g2,K,H_el,H_mag,H_gauss,H_gauss_closed_bot,H_gauss_closed_top,linkDim=max(20, ceil(Int,Nx*Ny*1.25)); rand_init=false, minsweeps=40, noise_sweeps=15, pass_init=false, init_state=0) # performs DMRG for given g². Returns energy and state
    if rand_init # option to start with a random initial state
        init = random_mps(sites;linkdims=linkDim)
    elseif pass_init # option to pass an initial state
        init = init_state
        linkDim = max(maxlinkdim(init),linkDim) # adjust initial link dimension based on link dimension of passed state
        linkDim = min(linkDim, 256) # however at most 256
    else # generates vaccuum state
        init = MPS(sites,["0" for n in 1:N])
    end
    
    obs = DMRGObserver(energy_tol=5E-5, minsweeps=minsweeps) # terminate DMRG when energy changes less than 5E-5 and at least minsweep (default is 50) sweeps have been performed
    maxdim = 10000 # limit maximum linkdim to avoid running out of ram
    energy,psi = dmrg([g2*H_el, 1/g2*H_mag,K*H_gauss,K*(H_gauss_closed_bot),K*(H_gauss_closed_top)], init ;nsweeps=10000,cutoff=1E-9, mindim=5, eigsolve_tol=1e-14,maxdim=Int.(reduce(vcat,[
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


function save_graph_links(state, name) # generates heatmap of expectation values of projector 1 - |0> <0|
    sites = []
    for y=1:Ny
        for x=1:Nx
            if x%2==0
                yoff=sqrt(3)/2
            else
                yoff=0
            end
            push!(sites,(3/2*(x-1),-(sqrt(3)*(y-1)+yoff)))
        end
    end
    p_sites = Point2f[sites...]

    right = (1.0,0.0)
    up = (-0.5,sqrt(3)/2)
    down = (-0.5,-sqrt(3)/2)

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

    linkwidth=50
    f = Figure()

    ax = Axis(f[1,1], autolimitaspect=1, width=3/2*Nx*linkwidth, height=sqrt(3)*(Ny+1/2)*linkwidth)
    hidedecorations!(ax)
    tightlimits!(ax)

    for i=1:N
        lines!(stack([p_sites[i],p_right[i]]), color=c_right[i],linewidth=10)
        lines!(stack([p_sites[i],p_up[i]]), color=c_up[i],linewidth=10)
        lines!(stack([p_sites[i],p_down[i]]), color=c_down[i],linewidth=10)
    end

    Colorbar(f[1,2], limits=(0,1), colormap=Reverse(:thermal))

    scatter!(p_sites, color="black", markersize=19, marker=:ltriangle)
    resize_to_layout!(f)

    save(name*".pdf", f, pdf_version="1.4")
end

function ground_sweep(g2s; carryover=false, w=1)
    l = length(g2s)
    # generate all hamiltonians and penalty terms
    H_el = get_H_el()
    H_mag = get_H_mag()
    H_gauss_closed_bot = get_H_gauss_closed_bottom()
    H_gauss_closed_top = get_H_gauss_closed_top()
    H_gauss = get_H_gauss()

    Ks = 100 .*[max(g2,1/g2) for g2 in g2s] # set constant for penalty terms based on parameter regime

    # initialize all the output lists
    gauss=[]
    closed_bottom=[]
    closed_top=[]
    variance=[]
    E_mags=[]
    E_els=[]
    num_rishons=[]

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
                state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed_bot,H_gauss_closed_top,noise_sweeps=15,rand_init=false)
                init = state
            else # any further iteration has minsweeps of 7 and takes previous state as initial state
                state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed_bot,H_gauss_closed_top,pass_init=true, init_state=init ,noise_sweeps=4, minsweeps=7)
                init = state
            end
        else # option for just independently using DMRG with random initial state
            state, energy = get_groundstate(g2s[i],Ks[i],H_el,H_mag,H_gauss,H_gauss_closed_bot,H_gauss_closed_top,noise_sweeps=25)
        end

        # calculate and save observables

        push!(gauss,inner(state',H_gauss,state))
        push!(closed_bottom,inner(state',H_gauss_closed_bot,state))
        push!(closed_top,inner(state',H_gauss_closed_top,state))
        push!(variance,inner(g2s[i]*H_el+1/g2s[i]*H_mag,state,g2s[i]*H_el+1/g2s[i]*H_mag,state)-energy^2)
        push!(Es,energy)
        push!(ex0s,expect(state,"p0"))
        push!(ex1s,expect(state,"p1"))
        push!(ex2s,expect(state,"p2"))
        push!(ex3s,expect(state,"p3"))

        push!(E_mags,inner(state',1/g2s[i]*H_mag,state))
        push!(E_els,inner(state',g2s[i]*H_el,state))
        push!(num_rishons,inner(state',H_el*4/6,state))

        save_graph_links(state, @sprintf("link_%ux%u_g2_%2.2f", Nx, Ny, g2s[i]))

        writedlm(@sprintf("AAA_g2_%ux%u%s", Nx, Ny, ".txt"), g2s)
        writedlm(@sprintf("AAA_E_%ux%u%s", Nx, Ny, ".txt"), Es)
        writedlm(@sprintf("AAA_gauss_%ux%u%s", Nx, Ny, ".txt"), gauss)
        writedlm(@sprintf("AAA_closed_bottom_%ux%u%s", Nx, Ny, ".txt"), closed_bottom)
        writedlm(@sprintf("AAA_closed_top_%ux%u%s", Nx, Ny, ".txt"), closed_top)
        writedlm(@sprintf("AAA_variance_%ux%u%s", Nx, Ny, ".txt"), variance)
        writedlm(@sprintf("AAA_ex0_%ux%u%s", Nx, Ny, ".txt"), ex0s)
        writedlm(@sprintf("AAA_ex1_%ux%u%s", Nx, Ny, ".txt"), ex1s)
        writedlm(@sprintf("AAA_ex2_%ux%u%s", Nx, Ny, ".txt"), ex2s)
        writedlm(@sprintf("AAA_ex3_%ux%u%s", Nx, Ny, ".txt"), ex3s)

        writedlm(@sprintf("AAA_num_rishon%ux%u%s", Nx, Ny, ".txt"), num_rishons)
        writedlm(@sprintf("AAA_E_mag%ux%u%s", Nx, Ny, ".txt"), E_mags)
        writedlm(@sprintf("AAA_E_el%ux%u%s", Nx, Ny, ".txt"), E_els)
    end
end


try
    ground_sweep(LinRange(8,0.5,31),carryover=true, w=parse(Int,ARGS[3])) # call julia dmrg_5xN.jl N_y w to start at the w-th g2 value. Usefull after crashes or timeouts
catch e
    ground_sweep(LinRange(8,0.5,31),carryover=true) # normal call that just runs through in order
end

#ground_sweep(LinRange(parse(Float64,ARGS[3]),parse(Float64,ARGS[3]),1),carryover=true) # call julia dmrg_5xN.jl
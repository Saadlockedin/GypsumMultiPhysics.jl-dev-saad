function Flux_Ficks_Array_Darcy!(f, u, edge, data)

    nk1 = edge.node[1]
    nl1 = edge.node[2]

     # Here calculate the multiComponent Diffusion 
     ϵ   = (data.E["ϵ"][nk1] + data.E["ϵ"][nl1])/2 #porosity
     τ   = (data.E["τ"][nk1] + data.E["τ"][nl1])/2 #tortuosity  
     #ΦPo = (data.E["pore_dia"][nk1] + data.E["pore_dia"][nl1])/2 #Pore dia
     #ΦPa = (data.E["part_dia"][nk1] + data.E["part_dia"][nl1])/2 #Part Dia
     #pm = Properties(ϵ,τ,ΦPo,ΦPa)
     # D_ij_e = (pm.ϵ/pm.τ) .* D_ij(data.E["tData"], T, P, data.E["tObj"].molwt)
     D_ije = D_ij(data.E["tData"], data.T, data.P, data.E["tObj"].molwt)
    for i = 1:data.E["nSp"]
        data.E["mF"][i] = (u[i,1]+u[i,2])/2.0;
    end

    totalMF = sum(data.E["mF"]);

    for i = 1:data.E["nSp"]
        data.E["mF"][i] = data.E["mF"][i]/totalMF;
    end

     # Knudsen Diffusion 
    # DK = (pm.ϵ/pm.τ) * D_Kn(data.E["tObj"].molwt, pm, T)
    # DK =  D_Kn(data.E["tObj"].molwt, pm, data.T)
    # Here we calculate the molecular diffusion 
    # Initialization
    suma = 0.0
    Dm = zeros(data.E["nSp"])
    Mm = 0.0
    alpha = zeros(data.E["nSp"])
    Deff = zeros(data.E["nSp"])

    # First loop: Calculate Molecular diffusion of each specie
    for i in 1:data.E["nSp"]
        suma = 0.0
        for j in 1:data.E["nSp"]
            if i != j
                suma += data.E["mF"][j]/ (D_ije[i,j])
            end
        end
        Dm[i] = (1.0 - data.E["mF"][i]) / (suma)
    end

    # Second loop: Calculate average molecular weight in the mixture 
    for i in 1:data.E["nSp"]
        Mm += data.E["mF"][i] * data.E["tObj"].molwt[i]
    end

    # Third loop: Calculate the effective diffusion of the specie ith in the mixture 
    for i in 1:data.E["nSp"]
        alpha[i] = 1.0 - sqrt(data.E["tObj"].molwt[i] / (Mm))
        Deff[i] = (pm.ϵ/pm.τ) * (1.0 / ((1.0 - alpha[i] * data.E["mF"][i]) / (Dm[i]) ))
    end


    μ = viscosity(data.E["tData"],data.T,data.E["tObj"].molwt,data.E["mF"])
    # Calculate ∇Pₜ
    Pₜ = 0.0; 
    for i = 1:data.E["nSp"]
        Pₜ = Pₜ + (data.R*data.T)*(u[i, 1] - u[i, 2])  
    end

    for i = 1:data.E["nSp"]
                
            #vh = (K/(μ*DK[i]))*Pₜ
            vh = (K/(μ))*Pₜ
            Bplus = fbernoulli(vh/Deff[i])
            Bminus = fbernoulli(-vh/Deff[i])
            f[i]= Deff[i]*((Bminus*u[i,1])-Bplus*u[i,2])
 
    end  

end

function Flux_Ficks_Array!(f, u, edge, data)

    nk1 = edge.node[1]
    nl1 = edge.node[2]

     # Here calculate the multiComponent Diffusion 
     ϵ   = (data.E["ϵ"][nk1] + data.E["ϵ"][nl1])/2
     τ   = (data.E["τ"][nk1] + data.E["τ"][nl1])/2
     ΦPo = (data.E["pore_dia"][nk1] + data.E["pore_dia"][nl1])/2
     ΦPa = (data.E["part_dia"][nk1] + data.E["part_dia"][nl1])/2
     pm = Properties(ϵ,τ,ΦPo,ΦPa)
     # D_ij_e = (pm.ϵ/pm.τ) .* D_ij(data.E["tData"], T, P, data.E["tObj"].molwt)
     D_ije = D_ij(data.E["tData"], data.T, data.P, data.E["tObj"].molwt)
    for i = 1:data.E["nSp"]
        data.E["mF"][i] = (u[i,1]+u[i,2])/2.0;
    end

    totalMF = sum(data.E["mF"]);

    for i = 1:data.E["nSp"]
        data.E["mF"][i] = data.E["mF"][i]/totalMF;
    end

     # Knudsen Diffusion 
    # DK = (pm.ϵ/pm.τ) * D_Kn(data.E["tObj"].molwt, pm, T)
     DK =  D_Kn(data.E["tObj"].molwt, pm, data.T)
    # Here we calculate the molecular diffusion 
    # Initialization
    suma = 0.0
    Dm = zeros(data.E["nSp"])
    Mm = 0.0
    alpha = zeros(data.E["nSp"])
    Deff = zeros(data.E["nSp"])

    # First loop: Calculate Molecular diffusion of each specie
    for i in 1:data.E["nSp"]
        suma = 0.0
        for j in 1:data.E["nSp"]
            if i != j
                suma += data.E["mF"][j]/ (D_ije[i,j])
            end
        end
        Dm[i] = (1.0 - data.E["mF"][i]) / (suma)
    end

    # Second loop: Calculate average molecular weight in the mixture 
    for i in 1:data.E["nSp"]
        Mm += data.E["mF"][i] * data.E["tObj"].molwt[i]
    end

    # Third loop: Calculate the effective diffusion of the specie ith in the mixture 
    for i in 1:data.E["nSp"]
        alpha[i] = 1.0 - sqrt(data.E["tObj"].molwt[i] / (Mm))
        Deff[i] = (pm.ϵ/pm.τ) * (1.0 / ((1.0 - alpha[i] * data.E["mF"][i]) / (Dm[i]) + 1.0 / (DK[i])))
    end


        
 # Calculate diffusive fluxes here 
    for i = 1:data.E["nSp"]
            f[i] = Deff[i]*(u[i,1]-u[i,2])
    end

end
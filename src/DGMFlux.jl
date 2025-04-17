function DGM_Diffusion(DK::Array{Float64,1}, Dij::Matrix{Float64}, Yᵢ)
    n = size(DK, 1)
    H = Matrix{Float64}(undef, n, n)
    for k=1:n
        sum = 0.0
        for l=1:n
            H[k,l] = -Yᵢ[k]/Dij[k,l]
            if k != l
                sum += Yᵢ[l]/Dij[k,l]
            end
        end
        H[k,k] = (1/DK[k]) + sum       
    end    
    H .= inv(H)
    return H
end

# This function will calculate the dusty gas diffusion 
# porosity :ϵ, tortuosity :τ, pore diamter and particle diameter should be given as an array
# Alternatively, use can use the Flux_DGM

function Flux_DGM_Array!(f, u, edge, E, R, T, P)
    #=
        E  = data.E; 
        R  = data.R
        T  = data.T;
    =#
     
        nk1 = edge.node[1]
        nl1 = edge.node[2]
    
         # Here calculate the multiComponent Diffusion 
         ϵ   = (E["ϵ"][nk1] + E["ϵ"][nl1])/2
         τ   = (E["τ"][nk1] + E["τ"][nl1])/2
         ΦPo = (E["pore_dia"][nk1] + E["pore_dia"][nl1])/2
         ΦPa = (E["part_dia"][nk1] + E["part_dia"][nl1])/2
    
         pm = Properties(ϵ,τ,ΦPo,ΦPa)
         D_ij_e = (pm.ϵ/pm.τ) .* D_ij(E["tData"], T, P, E["tObj"].molwt)
        for i = 1:E["nSp"]
            E["mF"][i] = (u[i,1]+u[i,2])/2.0;
        end
    
        totalMF = sum(E["mF"]);
    
        for i = 1:E["nSp"]
            E["mF"][i] = E["mF"][i]/totalMF;
        end
    
         # Knudsen Diffusion 
         DK = (pm.ϵ/pm.τ) * D_Kn(E["tObj"].molwt, pm, T)
         # Dusty Gas Diffusion 
         DGM = DGM_Diffusion(DK, D_ij_e, E["mF"])
         # Mixture viscosity
         μ = viscosity(E["tData"],T,E["tObj"].molwt,E["mF"]) 
        # end of diffusion calculations
    
        # Calculate ∇Pₜ
        Pₜ = 0.0; 
        for i = 1:E["nSp"]
            Pₜ = Pₜ + (R*T)*(u[i, 1] - u[i, 2])  
        end
            
     # Calculate diffusive and convective fluxes here 
        for i = 1:E["nSp"]
                
                sum_term = 0.0
                
                for j = 1:E["nSp"]
                   vh = (K/(μ*DK[i]))*Pₜ
                    Bplus = fbernoulli(vh)
                    Bminus = fbernoulli(-vh)
                    sum_term += DGM[i,j]*((Bminus*u[j,1])-Bplus*u[j,2])
                end
                f[i] = sum_term
        end
        return f
    end
    
import RootSolvers: NewtonsMethod, NewtonsMethodAD, SecantMethod, BrentsMethod

# KA.@print only accepts literal strings, so we must
# branch to print which method is being used.

#####
##### ρeq
#####
for rsm in (:NewtonsMethod, :NewtonsMethodAD, :SecantMethod, :BrentsMethod)
    @eval function print_warning_ρeq(::Type{<:$(rsm)}, args...)
        (ρ, e_int, q_tot, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(ρeq):\n",
            $(string(rsm)),
            "(ρ=$ρ, e_int=$e_int, q_tot=$q_tot, T=$T, maxiter=$maxiter, tol=$tol)",
            "\n"
        )
    end
    @eval function print_warning_hpq(::Type{<:$(rsm)}, args...)
        (h, p, q_tot, T_guess, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(hpq):\n",
            $(string(rsm)),
            "(h=$h, p=$p, q_tot=$q_tot, T_guess=$T_guess, T=$T, maxiter=$maxiter, tol=$tol) = ",
            "\n"
        )
    end
    @eval function print_warning_peq(::Type{<:$(rsm)}, args...)
        (p, e_int, q_tot, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(peq):\n",
            $(string(rsm)),
            "(p=$p, e_int=$e_int, q_tot=$q_tot, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
    @eval function print_warning_ρpq(::Type{<:$(rsm)}, args...)
        (ρ, p, q_tot, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(ρpq):\n",
            $(string(rsm)),
            "(ρ=$ρ, p=$p, q_tot=$q_tot, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
    @eval function print_warning_ρθq(::Type{<:$(rsm)}, args...)
        (ρ, θ_liq_ice, q_tot, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(ρθq):\n",
            $(string(rsm)),
            "(ρ=$ρ, θ_liq_ice=$θ_liq_ice, q_tot=$q_tot, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
    @eval function print_warning_pθq(::Type{<:$(rsm)}, args...)
        (p, θ_liq_ice, q_tot, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(pθq):\n",
            $(string(rsm)),
            "(p=$p, θ_liq_ice=$θ_liq_ice, q_tot=$q_tot, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
    @eval function print_warning_TᵥρRH(::Type{<:$(rsm)}, args...)
        (Tᵥ, RH, ρ, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(TᵥρRH):\n",
            $(string(rsm)),
            "(Tᵥ=$Tᵥ, RH=$RH, ρ=$ρ, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
    @eval function print_warning_ρθq_nonlinear(::Type{<:$(rsm)}, args...)
        (θ_liq_ice, ρ, q_tot, q_liq, q_ice, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(ρθq_nonlinear):\n",
            $(string(rsm)),
            "(θ_liq_ice=$θ_liq_ice, ρ=$ρ, q_tot=$q_tot, q_liq=$q_liq, q_ice=$q_ice, T=$T, maxiter=$maxiter, tol=$tol) =",
            "\n"
        )
    end
end

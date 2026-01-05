import RootSolvers: NewtonsMethod, NewtonsMethodAD, SecantMethod, BrentsMethod

# KA.@print only accepts literal strings, so we must
# branch to print which method is being used.

for rsm in (:NewtonsMethod, :NewtonsMethodAD, :SecantMethod, :BrentsMethod)
    for (IV, name, args) in (
        (ρeq, "ρeq", (:ρ, :e_int, :q_tot, :T, :maxiter, :tol)),
        (phq, "hpq", (:h, :p, :q_tot, :T_guess, :T, :maxiter, :tol)),
        (peq, "peq", (:p, :e_int, :q_tot, :T, :maxiter, :tol)),
        (pρq, "ρpq", (:ρ, :p, :q_tot, :T, :maxiter, :tol)),
        (ρθ_liq_ice_q, "ρθ_liq_ice_q", (:ρ, :θ_liq_ice, :q_tot, :T_guess, :T, :maxiter, :tol)),
        (pθ_liq_ice_q, "pθ_liq_ice_q", (:p, :θ_liq_ice, :q_tot, :T_guess, :T, :maxiter, :tol)),
    )
        @eval function print_warning(::Type{<:$(rsm)}, ::$IV, args...)
            ($(args...),) = args
            return KA.@print(
                $("maxiter reached in saturation_adjustment($name):\n"),
                $(string(rsm)),
                $(
                    let
                        parts = Any["("]
                        for (i, a) in enumerate(args)
                            push!(parts, "$(string(a))=")
                            push!(parts, a)
                            if i < length(args)
                                push!(parts, ", ")
                            end
                        end
                        push!(parts, ")" * (name == "ρeq" ? "" : " ="))
                        Expr(:string, parts...)
                    end
                ),
                "\n"
            )
        end
    end

    # Deprecated function, to be removed:
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

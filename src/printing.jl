import RootSolvers

# KA.@print only accepts literal strings, so we must
# branch to print which method is being used.

for rsm in (
    RootSolvers.NewtonsMethod,
    RootSolvers.NewtonsMethodAD,
    RootSolvers.SecantMethod,
    RootSolvers.BrentsMethod,
)
    for (IV, name, args) in (
        (ρe, "ρe", (:ρ, :e_int, :q_tot, :T, :maxiter, :tol)),
        (ph, "ph", (:h, :p, :q_tot, :T_guess, :T, :maxiter, :tol)),
        (pe, "pe", (:p, :e_int, :q_tot, :T, :maxiter, :tol)),
        (pρ, "pρ", (:ρ, :p, :q_tot, :T, :maxiter, :tol)),
        (
            ρθ_li,
            "ρθ_li",
            (:ρ, :θ_liq_ice, :q_tot, :T_guess, :T, :maxiter, :tol),
        ),
        (
            pθ_li,
            "pθ_li",
            (:p, :θ_liq_ice, :q_tot, :T_guess, :T, :maxiter, :tol),
        ),
    )
        @eval function print_warning(::Type{<:$(rsm)}, ::$IV, args...)
            ($(args...),) = args
            $(
                let
                    call_args = Any[]
                    push!(call_args, "maxiter reached in saturation_adjustment($name):\n")
                    push!(call_args, string(rsm))
                    push!(call_args, "(")
                    for (i, a) in enumerate(args)
                        push!(call_args, "$(string(a))=")
                        push!(call_args, a) # Symbol representing variable
                        if i < length(args)
                            push!(call_args, ", ")
                        end
                    end
                    push!(call_args, ")" * (name == "ρe" ? "" : " ="))
                    push!(call_args, "\n")
                    Expr(:macrocall, GlobalRef(KA, Symbol("@print")), LineNumberNode(@__LINE__), call_args...)
                end
            )
        end
    end

    # Deprecated function, to be removed:
    @eval function print_warning_ρθq_nonlinear(::Type{<:$(rsm)}, args...)
        (θ_liq_ice, ρ, q_tot, q_liq, q_ice, T, maxiter, tol) = args
        return KA.@print(
            "maxiter reached in saturation_adjustment(ρθq_nonlinear):\n",
            $(string(rsm)),
            "(θ_liq_ice=", θ_liq_ice, ", ρ=", ρ, ", q_tot=", q_tot, ", q_liq=", q_liq, ", q_ice=", q_ice, ", T=", T, ", maxiter=", maxiter, ", tol=", tol, ") =",
            "\n"
        )
    end
end

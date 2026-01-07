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
            return nothing
        end
    end

    # Deprecated function, to be removed:
    @eval function print_warning_ρθq_nonlinear(::Type{<:$(rsm)}, args...)
        return nothing
    end
end

module Utils

import Thermodynamics

"""
    _print(args...)

Print whatever is in `args`. Compatible with both CPU and GPU.

If CUDA is not loaded, then KernelAbstractions is not loaded which means
`Base.print` is called. If CUDA is loaded, then KernelAbstractions is loaded and
`KernelAbstractions.@print` is called.
"""
function _print(args...)
    if !isnothing(Base.get_extension(Thermodynamics, :KernelAbstractionsExt))
        ka_print(args...)
    else
        Base.print(args...)
    end
end

function ka_print end

end

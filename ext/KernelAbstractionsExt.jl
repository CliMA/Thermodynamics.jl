module KernelAbstractionsExt

import Thermodynamics.Utils as Utils
import KernelAbstractions as KA

Utils.ka_print(args...) = KA.@print(args...)

end

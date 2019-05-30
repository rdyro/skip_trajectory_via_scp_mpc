using Pkg, Convex

own_dir_name = dirname(@__FILE__)
convex_dir_name = dirname(pathof(Convex))

sol_path = joinpath(convex_dir_name, "solution.jl")
mdv_path = joinpath(convex_dir_name, "atoms/affine/multiply_divide.jl")

sol_patch = joinpath(own_dir_name, "patches", "solution.patch")
mdv_patch = joinpath(own_dir_name, "patches", "multiply_divide.patch")

sol_cmd = Cmd(["patch", "-f", "-r", "-", "--no-backup-if-mismatch",
               "--read-only=ignore", sol_path, sol_patch])
mdv_cmd = Cmd(["patch", "-f", "-r", "-", "--no-backup-if-mismatch",
               "--read-only=ignore", mdv_path, mdv_patch])

#sol_cmd = Cmd(["patch", "--help"])
#mdv_cmd = Cmd(["patch", "--help"])

try
  sol_run = run(sol_cmd)
catch
end
try
  mdv_run = run(mdv_cmd)
catch
end

sol_match = match(r"ROBERT DYRO", read(sol_path, String))
mdv_match = match(r"ROBERT DYRO", read(mdv_path, String))

if sol_match != nothing && mdv_match != nothing
  println("################################################")
  println("Patching finished successfully, re-import Convex")
  println("################################################")
else
  println("################################################")
  println("Didn't work, maybe it's already patched?")
  println("################################################")
end

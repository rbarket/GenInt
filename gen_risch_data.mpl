# Set the current directory to where this file is located if the read function below does not work
read "risch_funcs.mpl";  # loads all functions in the mpl file

N := 10000:  # number of pairs to generate
A := Array(1..N):

ext_args := [x, x+1, x^2, x^3, 1/x]:
extensions := [ln, exp]:

for extension in extensions do
    for arg in ext_args do
        for i from 1 to N do
            if i mod 1000 = 0 then print(i): end if:
            pair := TR_lin_gen(pick_int(), extension(arg), m=pick_int());
            A[i] := map(x -> convert(x, string), pair);
        end do:

# Create string based on extension and argument
    if arg = 1/x then
        str_name := cat(convert(extension, string), "(x^(-1))");
    else
        str_name := convert(extension(arg), string);
    end if:

# Export Data
    file_name := cat("generated_data_", str_name, ".json"):
    Export(file_name, A):
    end do:
end do:


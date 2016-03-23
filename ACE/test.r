Rcpp::sourceCpp('ACE_transform_cpponly.cpp')
system.time(ACE_transform_cpp())
print("OK")

# Rcpp::sourceCpp('ACE_std_modified.cpp')
# ACE_standard_cpp();
# print("OK")
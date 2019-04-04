# To regenerate the figures, you will need R, with the package circlize installed.
# For example, the circlize package can be installed as follows:
R --vanilla
install.packages(c("circlize"))
# You might want to use "sudo R --vanilla" if the package needs to be available for all users.
# On OSX, simply open a terminal, go in this directory, and type the following commands:


cd StrainComparisons  ; ./cmd.sh ; cd ..
cd HorizontalTransfer ; ./cmd.sh ; cd ..


# Note that during the work, internal names (T1, T2, TA, TB, TC1, TC2, TD, TE) were used for the templates.
# These were substituted by harmonized names (Experimental Conditions 1 to 8) in the final manuscript.
# The file  TemplateConversionList.xlsx  provides the conversion key.




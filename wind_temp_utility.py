import numpy as np

def give_code_coeffs(gamma,e0):
	TWF = (2./3.)*gamma*e0
	VWF = np.sqrt(2.*e0*(1-gamma))
	print "TWF: {}".format(TWF)
	print "VWF: {}".format(VWF)


give_code_coeffs(0.0,6.845)
# give_code_coeffs(0.25,6.845)
give_code_coeffs(0.5,6.845)
# give_code_coeffs(0.75,6.845)
# give_code_coeffs(0.95,6.845)

# give_code_coeffs(0.0,6.845)
# give_code_coeffs(0.1,6.845)
# give_code_coeffs(0.2,6.845)
# give_code_coeffs(0.3,6.845)
# give_code_coeffs(0.4,6.845)

# give_code_coeffs(0.0,6.845)
# give_code_coeffs(0.0,9.12666666667)
# give_code_coeffs(0.0,13.69)
# give_code_coeffs(0.5,6.845)
# give_code_coeffs(0.75,6.845)
# give_code_coeffs(0.95,6.845)


def give_thermal_coeff(gamma,VWF):
	TWF = gamma/(3-3*gamma) * VWF**2.
	print "TWF: {}".format(TWF)

# give_thermal_coeff(0.0,3.7)
# give_thermal_coeff(0.25,3.7)
# give_thermal_coeff(0.40,3.7)
# give_thermal_coeff(0.50,3.7)
# give_thermal_coeff(0.75,3.7)
# give_thermal_coeff(0.95,3.7)
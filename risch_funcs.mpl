# Perform Hermite reduction recursively until there is a resulting integral that has a square-free denominator
# expr: the expression to perform hermite reduction on
# theta_sub: the highest level field extension we are working in
# Output: the portion of the integral that is already integrated, and the remaining part of the expression to be integrated (that is square-free)
hermite := proc(expr, theta_sub)
	local P, Q, q_list, qk, k, T, t, diff_qk, diff_sigma, sub_tau, SF, cont, rat_part, log_part, remaining, sigma, tau;
	global x, theta;

	if numer(expr)=0 then  # Edge case when the log_part is 0
		return 0, 0;
	end if;
	
	# Step 1: Factor denom into SFF and define variables for diophantine eqn	
	Q := denom(expr);
	SF := eval(sqrfree(eval(Q, theta_sub=theta), theta), theta=theta_sub);  # Double eval because input the SF must be a polynomial in theta
	cont := SF[1];
	q_list := sort(SF[2], (x,y) -> x[2] < y[2]);  # sorts SF list based on the degree (we want k to be last)
	qk,k := q_list[-1][1], q_list[-1][2];
	
	# Check if denominator is SF
	if k=1 then		
		return 0, expr  # remaining=0, log_part=expr`
	end if;

	P := numer(expr)/cont;
	T := product(q_list[t][1]^(q_list[t][2]), t=1..numelems(q_list)-1);  # all factors from square-free except for qk

	# Step 2: Solve Diophantine eqn
	try  gcdex(eval(T*diff(qk,x), theta_sub=theta), 
		eval(qk, theta_sub=theta),
		eval(P, theta_sub=theta), 
		theta, 
		'sigma', 'tau');  # Solve diophantine equation for sigma, tau
	diff_sigma := diff(sigma,x);

	catch:
		print("no soln dio", expr);
	end try;
	sub_tau := subs(theta=theta_sub,tau);
	
	# Step 3: Plug into Hermite eqn
     rat_part := (sigma/(1-k))/(qk^(k-1));
     log_part := cont*(sub_tau + (diff_sigma/(k-1))*T)/simplify(Q/qk);
	
	remaining, log_part := hermite(log_part, theta_sub);  # Recursive call for when the log_part is not SF

	return eval(rat_part + remaining, theta=theta_sub), eval(log_part, theta=theta_sub);
	
end proc:

# Calculate the roots of the Trager-Rothstein Resultant given an expression and the top-level field extension
TR_roots := proc(expr, theta_sub)	
	local a, b, diff_b, Rz;
	global x, theta;

	# Assumption: denom is square-free	
	a := numer(expr);
	b := denom(expr);

	# Make denominator monic
	a /= lcoeff(b,theta_sub);
	b /= lcoeff(b,theta_sub);
	
	diff_b := eval(diff(b,x), theta_sub = theta);  # derivative b
	a,b := eval(a,theta_sub=theta), eval(b,theta_sub=theta);  # sub our extention with theta
	
	Rz := primpart(resultant(a - z*diff_b, b, theta), z); # resultant
	return {solve(Rz=0,z)}; #solving the resultant=0 is equivalent to the roots of the resultant
end proc:

# Use the Trager-Rothstein method to calculate the integral of an expression with a square-free denominator. 
# If the denominator is not SF, use the hermite reduction method first.
TR := proc(expr, theta_sub)
	local a, b, c, d, g, v, V, z, j, alpha, diff_b, Rz, fac_list, fac, fact, integral;
	global x, theta;

	# Assumption: denom is square-free	
	a := eval(numer(expr));
	b := eval(denom(expr));

	if degree(eval(a, theta_sub=theta), theta) >= degree(eval(b, theta_sub=theta), theta) then
		print("ERROR: degree of numerator is greater than degree of denominator", a, b);
		return FAIL;
	end if;

	# Make denominator monic
	a /= lcoeff(b,theta_sub);
	b /= lcoeff(b,theta_sub);
	diff_b := eval(diff(b,x), theta_sub = theta);  # derivative b
	a,b := eval(a,theta_sub=theta), eval(b,theta_sub=theta);  # sub our extention with theta
	
	Rz := primpart(resultant(a - z*diff_b, b, theta), z);
	fac_list := factors(Rz)[2];  # 1st part of list is the content which we don't need
	integral := 0;

	# Iteriate through each distinct factor (mult. doesn't matter)
	for fac in fac_list do  # fac[1] is the expr, fac[2] is the multiplicity (not needed)
		fact := fac[1];
		d := degree(fact, z);
		
		if d=0 then  # factor doesn't involve z
			next	
						
		elif d=1 then
			c := solve(fac[1]=0,z);

			
		#	if not type(c, integer) then  # root is not a constant => expr not elementary integrable				
		#		return FAIL;
		#	end if;

			v := gcdex(a-c*diff_b, b, theta);
			
			g := 0;
			if evalb(exp = op(0,theta_sub)) then  # add the extra term for exponential subcase
				g := -c*degree(v, theta)*op(1, theta_sub); # op is u in exp(u)
			end if;
			
			integral += c*log(v) + g;
			
		else
			alpha := RootOf(fact);
			V := gcd(a-alpha*diff_b, b);  # Tried using gcdex but has a hard time with RootOf numbers I believe.  
			V /= lcoeff(V);
			
			if d=2 then
				c := allvalues(alpha);  # same as solving fact=0 for z
				
				for j from 1 to 2 do
					v[j] := subs(alpha=c[j],V);
					integral += c[j]*log(v[j]);
				end do;
				
			end if;
				
		end if;
		
	end do:
	
	return eval(integral, theta=theta_sub);
end proc:

# Input is the set of roots from the TR algorithm
# For each root, set it equal to a constant and then isolate one of the symbolic coefficients
# Output is the condition on the symbolic coefficient to make the expression elementary integrable 
TR_soln := proc(TRroots)
local condn, i, ind, f, sol;
	condn := []:
	for i from 1 to numelems(TRroots) do
		ind := indets(TRroots[i], function): # all functions such as a(x) or a'(x)  

		for f in ind do
			if has(f,diff) or has(ind minus {f}, f) then # f is a derivative or f exists as a derivative in the indeterminates
				next;			
			fi:
			
			sol := solve(TRroots[i]=C[i], f); # the root is equal to a constant C
			condn := [op(condn),f=sol];
			break;
		od:
	od:

	return condn;
end proc:


# Generate a random polynomial in the previous field.
rand_coeff := proc(prev_field, min_exp, max_exp)
return randpoly(prev_field, expons=rand(min_exp..max_exp), coeffs=rand(-5..5))
end proc:

# Generate linear denominators in theta
lin_rand_denom := proc(num_terms)
	local term, terms, i;
	terms := [];
	
	for i from 1 to num_terms do
		
		term := 0;
		while not has(term, theta) do: # if coefficient of theta is 0, redo
			term := rand_coeff([x], rand(-1..0)(), rand(0..1)())*theta + rand_coeff([x], rand(-1..0)(), rand(0..1)()); # a(x)*theta + b(x)
		od:
		
		terms := [op(terms),term];  				
	end do:

	return terms;
end proc:

# Simple function for picking 1 with a 75% chance and picking 2 with a 25% chance, used for deciding num_terms and m in TR_lin_gen.
pick_int := proc()
    local r:
    r := rand(0...1.0)():
    if r <= 0.75 then return 1;
    else return 2;
    end if:
end proc:


# Generate an [integrand, integral] pair based on the TR-algorithm and Hermite reduction.
# num_terms: how many denominators to create for the partial fraction representation
# theta_sub: Which elementary extension we are working with (log(u) or exp(u), u in Q(x))
# m: the multiplicity of the denominators of each partial fraction. If m=1, then the denominator is square-free and only the TR-algorithm is used. If m>1, perform hermite reduction until the denominator is square-free and then use the TR-algorithm
TR_lin_gen := proc(num_terms, theta_sub, {m:=1})
local denoms, input, i, j, r, con, den, herm_part, trag_part, item, items, expr, integrand, integral;
integrand := 0;
while integrand=0 do
	# Create denominator (note that exponential extensions require theta does NOT divide all factors q_i
	denoms := lin_rand_denom(num_terms);	
	if has(theta_sub, exp) then  # check theta doesnt divide denominator for exp extensions
		do denoms := lin_rand_denom(num_terms); until not ormap(divide, denoms, theta);	
	end if;
	
	# Create partial fraction form given denominator
	input := 0:
	for i from 1 to numelems(denoms) do # iterate through the factors
		for j from 1 to m do 			
			input += a[i,j](x)/(denoms[i])^j;
		od;    		
	od:
	input := eval(input, theta=theta_sub);
	#Check if square-free (depends on parameter m)
	den := denom(input);

	herm_part := 0;
	if gcdex(eval(den,theta_sub=theta), eval(diff(den,x), theta_sub=theta), theta)<>1 then #if gcd(den, den')=1, expression is SF
		herm_part, trag_part := hermite(input, theta_sub);
	else
		trag_part := input;
	end if;
	
	# Solve for the symbolic coefficients
	r := TR_roots(trag_part, theta_sub);
	r := TR_soln(r);
	# Fill random values for the constants and free variables
	for con in r do
		
		items := indets(con) minus {x, lhs(con)};	# all symbolic coefficients and constants (a(x), a'(x), C[1], etc.). Subtract main variable and dependant coefficients	
		for item in items do
		
			if has(item, diff) then
				next; 
			elif type(item, indexed) then # is C_i
				r := eval(r, item=rand(-5..5)()/rand(1..5)()); # replace C_i with a random rational with max integer 5
			elif type(item, function) then # is a_i(x)
				expr := randpoly(x, coeffs=rand(-5..5), expons=rand(-1..1)); # replace symbolic coefficient with a polynomial in the prev. field
				r := eval(r, item=expr);
				r := [op(r), item=expr]; 	
			else:
				print("irrelevant",item);
			end if;S
		
		end do:
	end do:
	
	integrand := eval(input, r);
	herm_part, trag_part := hermite(integrand, theta_sub);
	integral := herm_part + TR(trag_part, theta_sub);
end do:

return [integrand, integral];
	 
end proc:
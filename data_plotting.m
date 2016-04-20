% This piece of code is just for the data plotting purpose 

% The elements in B preserve their columnwise ordering from A.
points = csvread('Extra_data_file.dat') ;
x_points = points(1,1);
y_points = points(1,2);
z_points = points(1,3);
free_stream_density  = points(1,4);
free_stream_velocity = points(1,5);

while true

		%plotting the residuals 
		Residuals =  csvread('Residual.dat');
		time = Residuals(:,1);
		density_res = Residuals(:,2);
		x_mom_res = Residuals(:,3);
		y_mom_res = Residuals(:,4);
		z_mom_res = Residuals(:,5);

		grid = csvread('Flat_plate.dat') ;
		param = csvread('conserved_variables.dat') ;
		x= reshape(grid(:,1),[y_points,x_points]) ;
		y= reshape(grid(:,2),[y_points,x_points]) ;
		density=  reshape(param(:,1),[y_points,x_points]);
		densityu= reshape(param(:,2),[y_points,x_points]);
		densityv= reshape(param(:,3),[y_points,x_points]);
		densityw= reshape(param(:,4),[y_points,x_points]);
		energy=   reshape(param(:,5),[y_points,x_points]);

		for i = 1:y_points
		    for j =	1:x_points	
		        u(i,j) = densityu(i,j)/density(i,j) ; 
		        v(i,j) = densityv(i,j)/density(i,j) ; 
		        w(i,j) = densityw(i,j)/density(i,j) ; 
		        velocity(i,j) = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)) ;
		        % mach(i,j) = velocity(i,j) / sqrt(1.4*287.14*temperature(i,j)) ;
				
				% non_dimensinalisation
				non_dim_x(i,j) = (x(i,j)-y(1,1))/(y(y_points,x_points)-y(1,1));
				non_dim_y(i,j) = (y(i,j)-y(1,1))/(y(y_points,x_points)-y(1,1));
				non_dim_den(i,j) = density(i,j)/free_stream_density;
				non_dim_vel(i,j) = velocity(i,j)/free_stream_velocity;
				% non_dim_press(i,j) = pressure(i,j)/pressure(y_points,x_points);
				% non_dim_temp(i,j) = temperature(i,j)/temperature(y_points,x_points) ;

		    end
		end

		%plotting the velocity at the exit
		figure(1)
		quiver(x,y,u,v) ;
		
		% plotting the residuals 
		pause(4) ;
		figure(2)
		plot(time,density_res,'LineWidth',1) ;
		title('Density residuals v/s time');

		
		% subplot(2,2,1,'replace')
		% plot(non_dim_den  (:,i), non_dim_y(:,i) ,'-o');
		% title('Non dim density vs y/L') ;

		% subplot(2,2,2,'replace')

  %       subplot(2,2,3,'replace')
  %       plot(non_dim_press(:,i), non_dim_y(:,i) ,'-go');
		% title('Non dim pressure vs y/L') ;

		% subplot(2,2,4,'replace')
		% plot(non_dim_temp (:,i), non_dim_y(:,i) ,'-bo') ;
		% title('Non dim temperature vs y/L')

		
		% subplot(,,5)
		% plot(mach(:,i), non_dim_y(:,i) ,'-o'); hold on ;
		% title('Non dim mach vs y/L')
		
		% h = surf(non_dim_x,non_dim_y) ;
		% plot (grid(:,1),grid(:,2),'-o') ;

		% figure2 
		% colormap hsv ;
		% colorbar ;
		% hcb=colorbar ;
		% shading interp;

		% subplot(2,2,1)
		% h = surf(non_dim_x,non_dim_y,mach) ;
		% t = title(hcb,'\bf {Mach}') ;
		
		% subplot(2,2,2)
		% h = surf(non_dim_x,non_dim_y,non_dim_vel) ;
		% t = title(hcb,'$\frac{V}{V_{\infty}}$') ;
		
		% subplot(2,2,3)
		% h = surf(non_dim_x,non_dim_y,non_dim_den) ;
		% t = title(hcb,'$\frac{\rho}{\rho_{\infty}}$') ;
		
		% subplot(2,2,4)
		% h = surf(non_dim_x,non_dim_y,non_dim_press) ;
		% t = title(hcb,'$\frac{P}{P_{\infty}}$') ;
		% % subplot(2,2,4)
		% % h = surf(non_dim_x,non_dim_y,non_dim_temp) ;
		% % t = title(hcb,'$\frac{T}{T_{\infty}}$') ;


		% t = title(hcb,'\bf {Mach}') ;
		% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf Mach');
		% title(' \bf Mach Number (M), flow over a bump')

		% t = title(hcb,'$\frac{\rho}{\rho_{\infty}}$') ;
		% set(t,'Interpreter','Latex');
		% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf \rho / \rho_{\infty}');
		% title(' \bf Nondimensional density(\rho/ \rho_{\infty}), flow over a bump')

		% t = title(hcb,'$\frac{V}{V_{\infty}}$') ;
		% set(t,'Interpreter','Latex');
		% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf V/V{\infty}');
		% title(' \bf Nondimensional velocity(V/V_{\infty}), flow over a bump')

		% title(' \bf Nondimensional pressure(P/P_{\infty}), flow over a bump')
		% t = title(hcb,'$\frac{P}{P_{\infty}}$') ;
		% set(t,'Interpreter','Latex');
		% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf P/P{\infty}');

		% t = title(hcb,'$\frac{T}{T_{\infty}}$') ;
		% set(t,'Interpreter','Latex');
		% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf T/T{\infty}');
		% title(' \bf Nondimensional temperature (T/T_{\infty}), flow over a bump')


	pause(4) ;
	close all ;
end

# SOS-solitary-waves
In this work optical solitary waves are generated through a variational approach using the symbiotyc organism search to optimize the parameters of the ansatz functions.

## Physical model: the generalized nonlinear Schrödinger equation

The NLSE is a very well known equation that has been studied in optics.  The canonical form of traditional one-dimensional NLSE used in classical field theory is shown below
$$  i\frac{\partial \Psi}{\partial t} + \frac{1}{2} \frac{\partial^2 \Psi}{\partial x^2} \pm \abs{\Psi}^2\Psi = 0 $$
where $i = \sqrt{-1}$ and $\Psi$ is the corresponding wave function.
It is worth noting that the term $\abs{\Psi}^2\Psi$ in equation \ref{NLSE} represents the nonlinear self phase modulation in an ideal focusing or defocusing Kerr medium with sign $+$ or $-$, respectively. For this case, the nonlinear is integrable in the sense that analytic solutions can be found and analyzed. However, most nonlinear physical systems of importance optics involve non-Kerr media such as saturation, photorefractive and competing media  \cite{Chen2012}. Therefore, more generally the NLSE equation can be stated in an an arbitrary nonlinear medium, that for the case of two dimensional scenario is
\begin{equation}
    i\frac{\partial \Psi}{\partial z} + \nabla^2 \Psi + F(r,|\Psi|^2)\Psi = 0  \,,
\end{equation}
where $F(r,|\Psi|^2)\Psi$ is the corresponding function that describe the effects of the nonlinear medium. In our physical model, $\Psi$ stands for the complex optical field, $r=(x,y)$ are the transverse spatial coordinates, and $z$ is the longitudinal (propagation) spatial coordinate, and $\nabla^2 \Psi=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}$ stands for the transverse Laplacian. Note that the generalized NLSE is written in dimensionless form.

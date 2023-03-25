#= 相对平衡位置做小振动的耦合双摆模型：
条件：
1、两相同单摆，摆长为a，摆锤质量m；两摆锤间以弹性系数k的轻质弹簧耦合，弹簧原长度与两摆悬点间距相等。
2、两摆均竖直时处于平衡位置，令此时的体系总势能V = 0。
3、体系自由度为2，取广义坐标为两摆各自与竖直方向的偏角θ、ϕ，对应广义动量: P₁ = m*a^2*θ、P₂ = m*a^2*ϕ
4、体系有2s(s=2)个简正频率，为 ω₁ = (g/a)^0.5, ω₂ = -(g/a)^0.5, ω₃ = (2k/m + g/a)^0.5, ω₄ = -(2k/m + g/a)^0.5
=#
using Plots, FFTW, DifferentialEquations

function CP!(du, u, Cnum, t)
    m, a, k, g = Cnum
    # Hamilton正则方程
    du[1] = u[2]/(m*a^2)
    du[2] = k*a^2*(u[3] - u[1]) - m*g*a*u[1]
    du[3] = u[4]/(m*a^2)
    du[4] = k*a^2*(u[1] - u[3]) - m*g*a*u[3]
end

function main()
    # 参数设置
    m = 1.0     
    a = 1.0
    k = 50.0
    g = 9.8
    dt = 0.001
    tmax = 20
    # 广义坐标、广义动量初值
    θ₀ = 0.04
    P₁₀ = 0.0
    ϕ₀ = 0.0
    P₂₀ = 0.0

    # 求解Hamilton正则方程
    tspan = (0, tmax)
    u₀ = [θ₀; P₁₀; ϕ₀; P₂₀]
    Cnum = [m a k g]
    prob = ODEProblem(CP!, u₀, tspan, Cnum)
    sol = solve(prob, Tsit5(), saveat=dt)  # saveat 指定t步长
    
    t = sol.t
    num_t = length(t)
    Msol = zeros(4, num_t)
    for i in 1:num_t
        Msol[:, i] = sol.u[i]
    end

    # 绘制相轨道
    # θ-P_θ
    fig1 = plot(Msol[1, :], Msol[2, :], 
        color=:red, 
        label="Pendulum_1",
        title=string("Phase spase", "(k=", k, ")"),
        xlabel="θ",
        ylabel="P_θ")
    # ϕ-P_ϕ
    fig2 = plot(Msol[3, :], Msol[4, :], 
        label="Pendulum_2",
        xlabel="ϕ ",
        ylabel="P_ϕ")
    # θ-ϕ
    fig3 = plot(Msol[1, :], Msol[3, :],
        color=:black,
        label=string("θ₀=", θ₀, " ϕ₀=", ϕ₀),
        xlabel="θ",
        ylabel="ϕ")
    fig4 = plot(fig1, fig2, fig3, 
        layout=(3, 1),
        size=(700,800))
    #savefig(fig4, "CoupledPendulum_phase.png")
end
main()

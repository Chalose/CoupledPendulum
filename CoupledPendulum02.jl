#= 耦合双摆模型：
条件：
1、两单摆，摆长为a，摆锤质量m₁,m₂ ；两摆锤间以弹性系数k的轻质弹簧耦合，弹簧原长度与两摆悬点间距相等。
2、两摆均竖直时处于平衡位置，令此时的体系总势能V = 0。
3、体系自由度为2，取广义坐标为两摆各自与竖直方向的偏角θ、ϕ，对应广义动量: P₁ = m₁*a^2*θ、P₂ = m₂*a^2*ϕ
4、小振动(m₁=m₂)时体系有2s(s=2)个简正频率，理论值为 ω₁ = (g/a)^0.5, ω₂ = -(g/a)^0.5, ω₃ = (2k/m + g/a)^0.5, ω₄ = -(2k/m + g/a)^0.5
=#
ENV["GKS_ENCODING"] = "utf-8"
using Plots, DifferentialEquations, FFTW
plotlyjs()
const m₁ = 1.0
const m₂ = 1.0
const a = 1.0
const k = 100.0
const g = 9.8


function CP!(du, u, p, t)
    # Hamilton正则方程
    du[1] = u[2]/(m₁*a^2)
    du[2] = -k*a^2*sin(u[1] - u[3]) - m₁*g*a*sin(u[1])
    du[3] = u[4]/(m₂*a^2)
    du[4] = k*a^2*sin(u[1] - u[3]) - m₂*g*a*sin(u[3])
end

function main()
    # 参数设置
    dt = 0.001     # 步长
    N = 100000     # 信号采样长度
    tmax = 100     # 终止时间(>= dt*N)

    # 广义坐标、广义动量初值
    θ₀ = π/2
    P₁₀ = 0.0
    ϕ₀ = -π/4
    P₂₀ = 0.0

    # 求解Hamilton正则方程
    tspan = (0, tmax)
    u₀ = [θ₀; P₁₀; ϕ₀; P₂₀]
    prob = ODEProblem(CP!, u₀, tspan)
    sol = solve(prob, Tsit5(), saveat=dt)  # saveat 指定t步长

    t = sol.t
    num_t = length(t)
    Msol = zeros(4, num_t)
    for i in 1:num_t
        Msol[:, i] = sol.u[i]
    end

    # θ, ϕ频谱分析
    f = 1/dt * Vector(0:Int(N/2))/N
    ω = 2π * f
    # 双边频谱
    tr1_θ = abs.(fft(Msol[1, 1:N]))/N
    tr1_ϕ = abs.(fft(Msol[3, 1:N]))/N
    # 单边频谱
    tr2_θ = tr1_θ[1:Int(N/2+1)]*2
    tr2_ϕ = tr1_ϕ[1:Int(N/2+1)]*2

    # 绘制相轨道, 频谱
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
    # 频谱
    fig5 = plot([ω ω], [tr2_θ tr2_ϕ],
        color=[:red :blue],
        xlabel="ω",
        ylabel="A",
        label=["θ" "ϕ"],
        linewidth=1.5)
    savefig(fig4, "CoupledPendulum02_phase.png")
    savefig(fig5, "CoupledPendulum02_omega.png")
end
main()

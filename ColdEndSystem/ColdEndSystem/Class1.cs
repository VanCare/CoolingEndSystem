//======================================================================
// 文件名:	ColdEndSystem

// 文件功能：电厂冷端系统（凝汽器，冷却塔）的计算：
//		固化输入设备、工质、环境等参数,
//		由热力计算程序提供排气焓值&排气流量， 
//		经过计算，输出凝汽器压力。
//		输出的凝汽器压力提供给热力计算程序进行迭代。
//		本DLL中的各个参数的单位均是国际标准单位。

//文件编写者：王凯(微信:wangkaiwangkaiwankai)

//文件修改时间：now
//======================================================================


using System;
using TurbinParameter;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

namespace ColdEndSystem
{

	/// <summary>
	/// 电厂冷端系统计算总程序
	/// ************************************************************************************
	/// 输入量：
	/// condenser	冷凝器数据，具体参看Condenser类的说明
	/// coolingTower	冷却塔数据，具体参看CoolingTower类的说明
	/// inTowerAir	大气数据，包括大气温度temperature、相对湿度relativeHumidity、压力pressure
	/// outTowerWater	出塔水（凝汽器入口水）数据：包括温度temperature（按环境湿球温度取值）、流量quantityOfFlow（由循环水泵计算，手动更改）、比热specificHeatCapacity
	/// inTurbineSteam	汽轮机排汽（凝汽器入口蒸汽）数据：包括排汽流量（来自热力计算程序）、排气焓（来自热力计算程序）
	/// outTurbineSteam	乏汽流出凝汽器部分数据，包含焓enthalpy
	/// ************************************************************************************
	/// 输出量：
	/// condenser.consenderPressure	凝汽器压力（提供给热力计算程序进行迭代计算）
	/// </summary>
	public class Program
	{
		public Condenser condenser;
		public CoolingTower coolingTower;
		public TurbineSteam inTurbineSteam;
		public TurbineSteam outTurbineSteam;
		public CoolingWater inTowerWater;
		public CoolingWater outTowerWater;
		public Air inTowerAir;
		public Air outTowerAir;

		public Program(Air inTowerAir, CoolingWater outTowerWater, TurbineSteam inTurbineSteam, Condenser condenser, CoolingTower coolingTower, TurbineSteam outTurbineSteam)
		{
			this.inTowerAir = inTowerAir;
			this.outTowerWater = outTowerWater;
			this.inTurbineSteam = inTurbineSteam;
			this.outTurbineSteam = outTurbineSteam;
			this.condenser = condenser;
			this.coolingTower = coolingTower;

			this.condenser.inTurbineSteam = this.inTurbineSteam;
			this.condenser.outTurbineSteam = this.outTurbineSteam;
			this.condenser.inCoolingWater = this.outTowerWater;
			this.condenser.outCoolingWater = this.inTowerWater;

			this.coolingTower.inCoolingWater = this.inTowerWater;
			this.coolingTower.outCoolingWater = this.outTowerWater;
			this.coolingTower.inAir = this.inTowerAir;
			this.coolingTower.outAir = this.outTowerAir;
		}

		public void CulationProcedure()
		{
			double i;
			do
			{
				i = condenser.consenderPressure;
				double j;
				do
				{
					j = coolingTower.outCoolingWater.temperature;
					condenser.Program();
					coolingTower.Program();
				}
				while (Math.Abs(j - coolingTower.outCoolingWater.temperature) > 0.1);   //
				condenser.outTurbineSteam.pressure = condenser.consenderPressure;
			}
			while (Math.Abs(i - condenser.consenderPressure) > 10000);	//
		}
	}

	/// <summary>
	/// 凝汽器数据
	/// ************************************************************************************
	/// 输入量：
	/// pipeDiameter	管内径
	/// pipeThickness	管壁厚
	/// pipeNumber	管输量
	/// pipeMeterial	管材
	/// flowPath	流程
	/// heatTransferArea	换热面积
	/// /// inCoolingWater	冷却水进入凝汽器部分数据，包含温度temperature（冷却塔计算）、流量quantityOfFlow（由循环水泵计算，手动更改）、比热specificHeatCapacity
	/// /// inTurbineSteam	乏汽进入凝汽器部分数据，包含流量quantityOfFlow、焓enthalpy
	/// /// outTurbineSteam	乏汽流出凝汽器部分数据，包含焓enthalpy
	/// ************************************************************************************
	/// 输出量：
	/// outCoolingWater	冷却水流出凝汽器数据，具体数据是温度temperature（提供给冷却塔）
	/// consenderPressure	凝汽器压力（提供给热力计算程序）
	/// ************************************************************************************
	/// 中间量：
	/// unjustedCoefficient	未修正凝汽器传热系数
	/// waterTemperatureCoefficient	水温修正系数
	/// pipeCoefficient	管材规格修正系数
	/// cleannessCoifficiente	清洁系数
	/// airLeakCoefficient	漏空气修正系数
	/// heatTransferCoefficient	修正凝汽器传热系数
	/// waterDeltaT	冷却水温升
	/// teminalTemperatureDiff	换热端差
	/// saturatedWaterTemperature	凝汽器饱和水温
	/// dispersionCoefficient	散质系数
	/// </summary>
	public class Condenser
	{
		public CoolingWater inCoolingWater;
		public CoolingWater outCoolingWater;
		public TurbineSteam inTurbineSteam;
		public TurbineSteam outTurbineSteam;

		public double pipeDiameter;
		public double pipeThickness;
		public int pipeNumber;
		public string pipeMeterial;
		public int flowPath;
		public double heatTransferArea;
		public double unjustedCoefficient;
		public double waterTemperatureCoefficient;
		public double pipeCoefficient;
		public double cleannessCoifficiente;
		public double airLeakCoefficient;
		public double heatTransferCoefficient;
		public double waterDeltaT;
		public double teminalTemperatureDiff;
		public double saturatedWaterTemperature;
		public double consenderPressure;

		public Condenser(double pipeDiameter,
						double pipeThickness,
						int pipeNumber,
						string pipeMeterial,
						int flowPath,
						double heatTransferArea)
		{
			this.pipeDiameter = pipeDiameter;
			this.pipeThickness = pipeThickness;
			this.pipeNumber = pipeNumber;
			this.pipeMeterial = pipeMeterial;
			this.flowPath = flowPath;
			this.heatTransferArea = heatTransferArea;
		}

		/// <summary>
		/// 冷却管内冷却水流速
		/// </summary>
		public void SolveWaterSpeed()
		{
			inCoolingWater.speedOfFlow = outCoolingWater.speedOfFlow = (4 * inCoolingWater.quantityOfFlow * flowPath) / (Math.PI * Math.Pow(pipeDiameter, 2.0) * pipeNumber);
		}

		/// <summary>
		/// 确定未修正传热系数
		/// </summary>
		public void SolveUnjustedCoefficient()
		{
			unjustedCoefficient = 1;
		}

		/// <summary>
		/// 水温修正系数
		/// </summary>
		public void SolveWaterTemperatureCoefficient()
		{
			waterTemperatureCoefficient = 1;
		}

		/// <summary>
		/// 管材规格修正系数
		/// </summary>
		public void SolvePipeCoefficient()
		{
			pipeCoefficient = 1;
		}

		/// <summary>
		/// 清洁系数
		/// </summary>
		public void SolveCleannessCoifficiente()
		{
			cleannessCoifficiente = 1;
		}

		/// <summary>
		/// 漏空气系数
		/// </summary>
		public void SolveAirLeakCoefficient()
		{
			airLeakCoefficient = 1;
		}

		/// <summary>
		/// 修正后的传热系数
		/// </summary>
		public void SolveHeatTransferCoefficient()
		{
			heatTransferCoefficient = (unjustedCoefficient * waterTemperatureCoefficient * pipeCoefficient * cleannessCoifficiente * airLeakCoefficient);
		}

		/// <summary>
		/// 循环水温升
		/// </summary>
		public void SolveWaterDeltaT()
		{
			waterDeltaT = (inTurbineSteam.quantityOfFlow * (inTurbineSteam.enthalpy - outTurbineSteam.enthalpy)) / (inCoolingWater.quantityOfFlow * inCoolingWater.specificHeatCapacity);
			outCoolingWater.temperature = inCoolingWater.temperature + waterDeltaT;
		}

		/// <summary>
		/// 换热端差
		/// </summary>
		public void SolveTeminalTemperatureDiff()
		{
			teminalTemperatureDiff = waterDeltaT / (Math.Exp((heatTransferCoefficient * heatTransferArea) / (inCoolingWater.specificHeatCapacity * inCoolingWater.quantityOfFlow)) - 1);
		}

		/// <summary>
		/// 凝汽器饱和水温
		/// </summary>
		public void SolveSaturatedWaterTemperature()
		{
			saturatedWaterTemperature = teminalTemperatureDiff + waterDeltaT + inCoolingWater.temperature;
		}

		/// <summary>
		/// 凝汽器压力
		/// </summary>
		public void SolveConsenderPressure()
		{
			consenderPressure = Water.T2P(saturatedWaterTemperature);	//
		}

		public void Program()
		{
			SolveWaterSpeed();
			SolveUnjustedCoefficient();
			SolveWaterTemperatureCoefficient();
			SolvePipeCoefficient();
			SolveCleannessCoifficiente();
			SolveAirLeakCoefficient();
			SolveHeatTransferCoefficient();
			SolveWaterDeltaT();
			SolveTeminalTemperatureDiff();
			SolveSaturatedWaterTemperature();
		}
	}

	/// <summary>
	/// 冷却塔数据
	/// ************************************************************************************
	/// 输入量：
	/// /// inCoolingWater	流入冷却塔的冷却水参数，需给出水温temperature，循环水流量quantityOfFlow,比热specificHeatCapacity
	/// /// inTurbineSteam	流入冷却塔的空气即大气的参数，需给出大气温度temperature、相对湿度relativeHumidity、压力pressure
	/// inletAreaRatio	冷却塔进风面积比
	/// packingResistanceCoefficient	填料阻力系数
	/// packingHeight	填料平均计算高度
	/// packingDepth	填料上缘至塔顶高度
	/// packingDiameter	填料面直径
	/// inletDiameter	进风口上缘直径
	/// inletHeight	进风口高度
	/// outletDiameter	冷却塔顶直径
	/// packingParamA	填料参数A
	/// packingParamM	填料参数M
	/// ************************************************************************************
	/// 输出量：
	/// outCoolingWater	流出冷却塔的冷却水参数，即水温temperature
	/// ************************************************************************************
	/// 中间量：
	/// airWaterRitio	空气和水的质量比
	/// midAir	填料断面的空气参数
	/// averageAir 温度为进塔空气、出塔空气平均温度时的空气参数
	/// saturatedHighEnthalpy 空气温度为进塔水温时饱和空气的比焓
	/// saturatedAverageEnthalpy	空气温度为平均水温时饱和空气的比焓
	/// saturatedLowEnthalpy	空气温度为出塔水温时饱和空气的比焓
	/// airResistance	通风阻力
	/// Omega	特征数Omega
	/// coolingNumberN	冷却数N
	/// </summary>
	public class CoolingTower
	{

		public CoolingWater inCoolingWater;
		public CoolingWater outCoolingWater;
		public Air inAir;
		public Air outAir;
		public Air midAir;
		public Air averageAir;

		public double inletAreaRatio;
		public double packingResistanceCoefficient;
		public double packingHeight;
		public double packingDepth;
		public double packingDiameter;
		public double inletDiameter;
		public double inletHeight;
		public double outletDiameter;
		public double packingParamA;
		public double packingParamM;

		public double airWaterRitio;
		public double saturatedHighEnthalpy;
		public double saturatedAverageEnthalpy;
		public double saturatedLowEnthalpy;
		public double dispersionCoefficient;
		public double resistanceCoefficientAll;
		public double resistanceCoefficientA;
		public double resistanceCoefficientB;
		public double resistanceCoefficientE;
		public double airResistance;
		public double Omega;
		public double coolingNumberN;

		public CoolingTower(double inletAreaRatio,
							double packingResistanceCoefficient,
							double packingHeight,
							double packingDepth,
							double packingDiameter,
							double inletDiameter,
							double inletHeight,
							double outletDiameter,
							double packingParamA,
							double packingParamM)
		{
			this.inletAreaRatio = inletAreaRatio;
			this.packingResistanceCoefficient = packingResistanceCoefficient;
			this.packingHeight = packingHeight;
			this.packingDepth = packingDepth;
			this.packingDiameter = packingDiameter;
			this.inletDiameter = inletDiameter;
			this.inletHeight = inletHeight;
			this.outletDiameter = outletDiameter;
			this.packingParamA = packingParamA;
			this.packingParamM = packingParamM;
		}
		
		/// <summary>
		/// 求冷却塔出口水温计算值
		/// </summary>
		public void Program()
		{
			SolveInAir();
			SolveSaturatedHighEnthalpy();
			SolveResistanceCoefficientA();
			SolveResistanceCoefficientE();
			RenewMidAirSpeed();
			do
			{
				SolveQuantityOfAir();
				SolveAirWaterRitio();
				RenewOutCoolingWaterTemperature();
				SolveEnthalpy();
				SolvePreciseAir();
				SolveCoolingNumberN();
				SolveCharacteristicNumberOmega();
			}
			while (JudgeAbs());
		}

		/// <summary>
		/// 根据 大气温度湿度压力 计算 进入空气的焓值和密度
		/// </summary>
		public void SolveInAir()
		{
			TemperatureHumidityPressure2Enthalpy(inAir);
			TemperatureHumidityPressure2Density(inAir);
		}

		/// <summary>
		/// 根据 进塔水温和大气压力 计算 温度为进塔水温时饱和空气的比焓
		/// </summary>
		public void SolveSaturatedHighEnthalpy()
		{
			saturatedHighEnthalpy = 1;  //
		}

		/// <summary>
		/// 更新填料断面的气流速度
		/// </summary>
		public void RenewMidAirSpeed()
		{
			if (midAir.speedOfFlow == 0)
			{
				midAir.speedOfFlow = 1.0;
			}
			else
			{
				midAir.speedOfFlow = Math.Pow(((4 * 9.8 * (packingDepth + 0.5 * packingHeight) * (inAir.density - outAir.density)) / (resistanceCoefficientAll * (inAir.density + outAir.density))), 0.5);
			}
		}

		/// <summary>
		/// 根据填料断面气流速度计算空气流量
		/// </summary>
		public void SolveQuantityOfAir()
		{
			inAir.quantityOfFlow = midAir.quantityOfFlow = outAir.quantityOfFlow = Math.PI * packingDiameter * packingDiameter * inAir.density * midAir.speedOfFlow / 4;
		}

		/// <summary>
		/// 根据空气流量计算气水比
		/// </summary>
		public void SolveAirWaterRitio()
		{
			airWaterRitio = inAir.quantityOfFlow / inCoolingWater.quantityOfFlow;
		}

		/// <summary>
		/// 更新出塔水温
		/// </summary>
		public void RenewOutCoolingWaterTemperature()
		{
			if (outCoolingWater.temperature == 0)
			{
				outCoolingWater.temperature = inCoolingWater.temperature - 8;
				do
				{
					outCoolingWater.temperature += 1;
					SolveEnthalpy();
				}
				while (JudgeDeltaEnthalpy());
			}
			else
			{
				if (coolingNumberN > Omega)
				{
					outCoolingWater.temperature += 0.5 * (inCoolingWater.temperature - outCoolingWater.temperature);
				}
				else
				{
					outCoolingWater.temperature -= 0.5 * (inCoolingWater.temperature - outCoolingWater.temperature);
				}
			}
		}

		/// <summary>
		/// 计算冷却数N对应的剩余四个焓值：outAir.enthalpy、averageAir.enthalpy、saturatedAverageEnthalpy、saturatedLowEnthalpy
		/// </summary>
		public void SolveEnthalpy()
		{
			Temperature2LatentHeatOfVaporization(outCoolingWater);
			dispersionCoefficient = 1 - inCoolingWater.specificHeatCapacity * outCoolingWater.temperature / outCoolingWater.latentHeatOfVaporization;
			outAir.enthalpy = inAir.enthalpy + inCoolingWater.specificHeatCapacity * (inCoolingWater.temperature - outCoolingWater.temperature) / (dispersionCoefficient * airWaterRitio);
			averageAir.enthalpy = (inAir.enthalpy + outAir.enthalpy) / 2;
			saturatedAverageEnthalpy = 1;   //
			saturatedLowEnthalpy = 1;   //
		}

		/// <summary>
		/// 迭代计算总阻力系数、填料断面气流速度等参数
		/// </summary>
		public void SolvePreciseAir()
		{
			double i;
			do
			{
				i = midAir.speedOfFlow;
				RenewPressureAndDensity();
				SolveResistanceCoefficientAll();
				RenewMidAirSpeed();
				SolveAirResistance();
			}
			while (Math.Abs(i - midAir.speedOfFlow) > 0.05);    //
		}

		/// <summary>
		/// 计算冷却数N
		/// </summary>
		public void SolveCoolingNumberN()
		{
			coolingNumberN = inCoolingWater.specificHeatCapacity * (inCoolingWater.temperature - outCoolingWater.temperature) / 6 * (1 / (saturatedLowEnthalpy - inAir.enthalpy) + 4 / (saturatedAverageEnthalpy - inAir.enthalpy) + 1 / (saturatedHighEnthalpy - outAir.enthalpy));
		}

		/// <summary>
		/// 计算特性数Omega
		/// </summary>
		public void SolveCharacteristicNumberOmega()
		{
			Omega = packingParamA * Math.Pow(airWaterRitio, packingParamM);
		}

		/// <summary>
		/// 判断N和Omega之差的绝对值是否大于0.01，大于返回true继续迭代，小于返回false推出迭代。
		/// </summary>
		/// <returns></returns>
		public bool JudgeAbs()
		{
			if (Math.Abs(coolingNumberN - Omega) > 0.01)
			{
				return true;
			}
			else
			{
				return false;
			}
		}


		/// <summary>
		/// 根据空气温度湿度压力求焓值
		/// </summary>
		/// <returns>空气焓值</returns>
		public void TemperatureHumidityPressure2Enthalpy(Air air)
		{
			air.enthalpy = 1;	//
		}

		/// <summary>
		/// 根据空气温度湿度压力求密度
		/// </summary>
		/// <returns>空气密度</returns>
		public void TemperatureHumidityPressure2Density(Air air)
		{
			air.density = 1;	//
		}

		/// <summary>
		/// 根据水温计算汽化潜热
		/// </summary>
		/// <returns>汽化潜热</returns>
		public void Temperature2LatentHeatOfVaporization(CoolingWater water)
		{
			water.latentHeatOfVaporization = 1;	//
		}

		/// <summary>
		/// 修正更新出塔空气的压力，同时校正密度
		/// </summary>
		public void RenewPressureAndDensity()
		{
			if (outAir.pressure == 0)
			{
				outAir.pressure = inAir.pressure = 1;   //
			}
			else
			{
				outAir.pressure = 1;    //划掉（修正h2对应的压力(第五页)，已知相对湿度为1温度焓）
			}
			outAir.density = 1; //
			averageAir.density = 1; //
		}

		/// <summary>
		/// 判断冷却数N公式中的三个焓差都大于0，防止迭代不收敛
		/// </summary>
		/// <returns></returns>
		public bool JudgeDeltaEnthalpy()
		{
			if (saturatedLowEnthalpy - inAir.enthalpy > 0 &&
				saturatedAverageEnthalpy - inAir.enthalpy > 0 &&
				saturatedHighEnthalpy - outAir.enthalpy > 0)
			{
				return false;
			}
			else
			{
				return true;
			}
		}

		/// <summary>
		/// 计算总阻力系数
		/// </summary>
		public void SolveResistanceCoefficientAll()
		{
			SolveResistanceCoefficientB();
			resistanceCoefficientAll = resistanceCoefficientA + resistanceCoefficientB + resistanceCoefficientE;
		}
		/// <summary>
		/// 计算阻力系数A
		/// </summary>
		public void SolveResistanceCoefficientA()
		{
			resistanceCoefficientA = (1 - 3.47 * inletAreaRatio + 3.65 * Math.Pow(inletAreaRatio,2)) * (85+2.51*packingResistanceCoefficient-0.206*Math.Pow(packingResistanceCoefficient,2)+0.00962*Math.Pow(packingResistanceCoefficient,3));
		}

		/// <summary>
		/// 计算阻力系数B
		/// </summary>
		public void SolveResistanceCoefficientB()
		{
			resistanceCoefficientB = 6.72 + 0.654 * packingDiameter + 3.5 * (4 * inCoolingWater.quantityOfFlow / (Math.PI * Math.Pow(packingDiameter, 2))) + 1.43 * midAir.speedOfFlow - 60.61 * inletAreaRatio - 0.36 * midAir.speedOfFlow * packingDiameter;
		}

		/// <summary>
		/// 计算阻力系数E
		/// </summary>
		public void SolveResistanceCoefficientE()
		{
			resistanceCoefficientE = Math.Pow((packingDiameter/outletDiameter), 4);
		}

		/// <summary>
		/// 计算冷却塔通风阻力
		/// </summary>
		public void SolveAirResistance()
		{
			airResistance = resistanceCoefficientAll * (inAir.density + outAir.density) * Math.Pow(midAir.speedOfFlow, 2) / 4;
		}
	}

	/// <summary>
	/// 工质抽象类，其子类包括蒸汽，冷却水，空气。  抽象特征为：拥有压力，温度，焓值，流量，比热等字段。
	/// </summary>
	public abstract class WorkingMedium
	{
		public double temperature;  //temperature
		public double pressure;
		public double quantityOfFlow;   //流量
		public double enthalpy;  //焓
		public double specificHeatCapacity; //比热容
		public double speedOfFlow;  //流速
		public double density;
		public WorkingMedium(double temperature = 0, double pressure = 0, double quantityOfFlow = 0, double enthalpy = 0, double specificHeatCapacity = 0, double speedOfFlow = 0, double density = 0)
		{
			this.temperature = temperature;
			this.pressure = pressure;
			this.quantityOfFlow = quantityOfFlow;
			this.enthalpy = enthalpy;
			this.specificHeatCapacity = specificHeatCapacity;
			this.speedOfFlow = speedOfFlow;
			this.density = density;
		}
	}

	/// <summary>
	/// 冷却水数据
	/// </summary>
	public class CoolingWater : WorkingMedium
	{
		public double latentHeatOfVaporization;
		public CoolingWater(double temperature, double pressure, double quantityOfFlow, double enthalpy, double specificHeatCapacity, double latentHeatOfVaporization = 0) : base(temperature, pressure, quantityOfFlow, enthalpy, specificHeatCapacity)
		{
			this.latentHeatOfVaporization = latentHeatOfVaporization;
		}
	}

	/// <summary>
	/// 透平蒸汽数据
	/// </summary>
	public class TurbineSteam : WorkingMedium
	{
		public TurbineSteam(double temperature, double pressure, double quantityOfFlow, double enthalpy, double specificHeatCapacity) : base(temperature, pressure, quantityOfFlow, enthalpy, specificHeatCapacity)
		{

		}
	}

	/// <summary>
	/// 空气数据
	/// </summary>
	public class Air : WorkingMedium
	{
		public double relativeHumidity;

		public Air(double temperature, double pressure, double quantityOfFlow, double enthalpy, double specificHeatCapacity, double density, double relativeHumidity = 0) : base(temperature, pressure, quantityOfFlow, enthalpy, specificHeatCapacity, density)
		{
			this.relativeHumidity = relativeHumidity;
		}
	}
}
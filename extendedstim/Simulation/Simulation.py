"""
Module: Simulation
"""
import multiprocessing
import os
import time
from multiprocessing import Pool
import numpy as np
import pymatching
import stim
import stimbposd
from extendedstim.Circuit.Circuit import Circuit


##  APPENDIX: Try to import tesseract decoder
try:
    from tesseract_decoder import tesseract
    import tesseract_decoder
except ImportError:
    print("Tesseract decoder not found")
    class tesseract:
        def __init__(self):
            pass
        @staticmethod
        def TesseractConfig(dem, pqlimit, det_orders, det_beam=50, beam_climbing=False, no_revisit_dets=True):
            return Config(dem, det_beam)
    class Config:
        def __init__(self, dem: stim.DetectorErrorModel, det_beam):
            self.dem = dem
            self.det_beam = det_beam
        def compile_decoder(self):
            return stimbposd.bp_osd.BPOSD(model=self.dem, bp_method='min_sum', max_bp_iters=self.det_beam)


sample_number_now = multiprocessing.Value('i', 0)
error_number_now = multiprocessing.Value('i', 0)


##  CHAPTER: ====Simulation====
def Simulation(circuit: Circuit, p_noise_list, p_measure_list, allowed_variant_ratios, decoder_name, decoder_params=None):
    """""        
    Tesseract decoder参数说明：
    {
    pqlimit：int对象，次数
    det_beam：int对象，分支查询数目
    beam_climbing：bool对象，是否梯度方法
    num_det_orders=num_det_orders，不知道什么阶数
    no_revisit_dets=no_revisit_dets，是否避免重复访问
    }
    ----------------------------------
    BP-OSD decoder参数说明：
    {
    bp_method：str对象，BP更新方法
    1. "product_sum"
    2. "min_sum"
    3. "min_sum_log"
    max_bp_iters：int对象，BP迭代次数，by default，50
    osd_method：str对象
        1. "osd_0"
        2. "osd_e": exhaustive OSD
        3. "osd_cs": combination-sweep OSD
        4. "osd_e"
    osd_order：int对象，OSD置换数目
    }
    """""

    prototype = circuit.prototype
    group_number = os.cpu_count() // 2

    logical_error_rates = []
    variants = []
    for i in range(len(p_noise_list)):
        dem = prototype2dem(prototype, p_noise_list[i], p_measure_list[i])
        start_time = time.time()
        task(dem, decoder_name, 20, decoder_params=decoder_params, test=True)
        end_time = time.time()
        duration = (end_time - start_time) / 100
        print('解码一个样本的时间为：', duration, '秒')
        if 50 / duration > 500000:
            sample_in_one_group_number = 500000
        elif 50 / duration < 10:
            sample_in_one_group_number = 10
        else:
            sample_in_one_group_number = int(50 / duration)

        lock = multiprocessing.Manager().Lock()
        global sample_number_now
        global error_number_now
        sample_number_now.value = 0
        error_number_now.value = 0
        with Pool(processes=group_number) as pool:
            results = [pool.apply_async(task, args=(dem, decoder_name, sample_in_one_group_number, lock, allowed_variant_ratios[i], decoder_params)) for _
                       in range(group_number)]
            _ = [result.get() for result in results]
        error_rate = error_number_now.value / sample_number_now.value
        sigma = np.sqrt((error_number_now.value / sample_number_now.value - (
                error_number_now.value / sample_number_now.value) ** 2) / sample_number_now.value) * 2 * 3.2905
        logical_error_rates.append(error_rate)
        variants.append(sigma)

    return p_noise_list,p_measure_list,logical_error_rates, variants


##  CHAPTER: ====Prototype to DEM====
def prototype2dem(prototype, p_noise, p_measure):
    dem_str = ''
    dem_temp = []
    for i, temp in enumerate(prototype):
        if temp[2] == 'DEP1':
            p_temp = (1 - np.sqrt(1 - 4 * p_noise / 3)) / 2
        elif temp[2] == 'DEP2':
            p_temp = (1 / 2) * (1 - (1 - 16 * p_noise / 15) ** (1 / 8))
        elif temp[2] == 'ONE_ERROR':
            p_temp = p_noise
        elif temp[2] == 'M_ERROR':
            p_temp = p_measure
        else:
            raise ValueError
        flag = True
        for j, item in enumerate(dem_temp):
            if temp[1] == item[1]:
                item[0] = p_temp * (1 - item[0]) + (1 - p_temp) * item[0]
                flag = False
                break
        if flag:
            dem_temp.append([p_temp, temp[1]])

    for temp in dem_temp:
        dem_str += ('\n' + f'error({temp[0]})' + temp[1])
    return stim.DetectorErrorModel(dem_str)


##  CHAPTER: ====Task====
def task(dem, decoder_name, number, lock=None, error_bar=None, decoder_params=None, test=False):
    sampler = dem.compile_sampler()  # 采样器

    if decoder_name == 'bposd':
        max_bp_iters = 50
        bp_method = 'min_sum'
        osd_method = 'osd_e'
        osd_order = 60
        if decoder_params is not None:
            if 'max_bp_iters' in decoder_params:
                max_bp_iters = decoder_params['max_bp_iters']
            if 'bp_method' in decoder_params:
                bp_method = decoder_params['bp_method']
            if 'osd_method' in decoder_params:
                osd_method = decoder_params['osd_method']
            if 'osd_order' in decoder_params:
                osd_order = decoder_params['osd_order']
        decoder = stimbposd.bp_osd.BPOSD(model=dem, bp_method=bp_method, max_bp_iters=max_bp_iters, osd_method=osd_method, osd_order=osd_order)

    elif decoder_name == 'tesseract':
        det_beam = 15
        beam_climbing = True
        no_revisit_dets = True
        num_det_orders = 16

        if decoder_params is not None:
            if 'det_beam' in decoder_params:
                det_beam = decoder_params['det_beam']
            if 'beam_climbing' in decoder_params:
                beam_climbing = decoder_params['beam_climbing']
            if 'num_det_orders' in decoder_params:
                num_det_orders = decoder_params['num_det_orders']
            if 'no_revisit_dets' in decoder_params:
                no_revisit_dets = decoder_params['no_revisit_dets']

        config = tesseract.TesseractConfig(dem=dem, pqlimit=200_000, det_beam=det_beam, beam_climbing=beam_climbing,
                                           det_orders=tesseract_decoder.utils.build_det_orders(
                                               dem=dem,
                                               num_det_orders=num_det_orders,
                                               method=tesseract_decoder.utils.DetOrder.DetIndex,
                                           ),
                                           no_revisit_dets=no_revisit_dets, )
        decoder = config.compile_decoder()

    elif decoder_name == 'matching':
        decoder = pymatching.Matching()
        decoder.from_detector_error_model(dem)
    else:
        raise NotImplementedError

    if test:
        detector_data, obs_data, error_data = sampler.sample(shots=number)  # 样本
        valid_detector_data = detector_data[
            np.where(np.logical_or(np.any(obs_data, axis=1), np.any(detector_data, axis=1)))[0]]  # 有效样本
        valid_obs_data = obs_data[np.where(np.logical_or(np.any(obs_data, axis=1), np.any(detector_data, axis=1)))[0]]  # 有效样本
        predictions = decoder.decode_batch(valid_detector_data)  # 解码器对每个样本的预测

        num_errors = 0
        for shot in range(len(valid_detector_data)):
            actual_for_shot = valid_obs_data[shot]
            predicted_for_shot = predictions[shot]
            if not np.array_equal(actual_for_shot, predicted_for_shot):
                num_errors += 1

    else:
        while True:
            global sample_number_now
            global error_number_now
            detector_data, obs_data, error_data = sampler.sample(shots=number)  # 样本
            valid_detector_data = detector_data[
                np.where(np.logical_or(np.any(obs_data, axis=1), np.any(detector_data, axis=1)))[0]]  # 有效样本
            valid_obs_data = obs_data[np.where(np.logical_or(np.any(obs_data, axis=1), np.any(detector_data, axis=1)))[0]]  # 有效样本
            predictions = decoder.decode_batch(valid_detector_data)  # 解码器对每个样本的预测

            num_errors = 0
            for shot in range(len(valid_detector_data)):
                actual_for_shot = valid_obs_data[shot]
                predicted_for_shot = predictions[shot]
                if not np.array_equal(actual_for_shot, predicted_for_shot):
                    num_errors += 1
            lock.acquire()
            sample_number_now.value += number
            error_number_now.value += num_errors
            error_rate = error_number_now.value / sample_number_now.value
            sigma = np.sqrt((error_number_now.value / sample_number_now.value - (
                    error_number_now.value / sample_number_now.value) ** 2) / sample_number_now.value)
            number_temp = sample_number_now.value
            lock.release()
            print(f'number: {number_temp}, error_rate: {error_rate:.6f}, sigma: {sigma:.6f}')
            if sigma < error_bar * error_rate:
                break

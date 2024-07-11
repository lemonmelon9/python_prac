class topography:
    
    '''
    Args:
        filepath : str
            Name of the Nanonis spectrum file to be loaded.
    
    Attributes (name : type):
        fname : str
            The name of the file excluding its containing directory.
        header : dict
            Header information of spectrum data.
        signals : dict
            Measured values in spectrum data.

    Methods:
        get_z(self, processing = 'raw', scan_direction = 'fwd')
            Parameters:
            processing : str
                The image processing.
                Possible parameters is the following: 
                    'raw', 'subtract average', 'subtract linear fit', 'subtract parabolic fit', 'differentiate'
            scan_direction : str
                The direction of scan.
                Possible parameters is the following: 
                    'fwd', 'bwd'
            
            Returns the two-dimensional array
            that represents the topographic data
            specifically processed for scanning in the direction of scan_direction.
            
            The detailed processing is carried out through the following five methods:
                raw, subtract_average, subtract_linear_fit, subtract_parabolic_fit, differentiate.
            
        raw (self, scan_direction)
            Returns the two-dimensional array
            containing the raw z data.
        
        subtract_average (self, scan_direction)
            Returns the two-dimensional array
            containing the z data processed through subtract average.
        
        subtract_linear_fit (self, scan_direction)
            Returns the two-dimensional array
            containing the z data processed through subtract linear fit.
        
        subtract_parabolic_fit (self, scan_direction)
            Returns the two-dimensional array
            containing the z data processed through subtract parabolic fit.
        
        differentiate (self, scan_direction)
            Returns the two-dimensional array
            containing the dz/dx data, in which x represents the fast scan axis.
    '''
    
    def __init__(self, filepath):
        import nanonispy as nap
        self.fname = nap.read.Scan(filepath).fname.split('/')[-1]
        self.header = nap.read.Scan(filepath).header
        self.signals = nap.read.Scan(filepath).signals
    
    def get_z(self, processing = 'raw', scan_direction = 'fwd'): # 'fwd' or 'bwd'
        if processing == 'raw':
            return self.raw(scan_direction)
        elif processing == 'subtract average':
            return self.subtract_average(scan_direction)
        elif processing == 'subtract linear fit':
            return self.subtract_linear_fit(scan_direction)
        elif processing == 'subtract linear fit xy':
            return self.subtract_linear_fit_xy(scan_direction)
        elif processing == 'subtract parabolic fit':
            return self.subtract_parabolic_fit(scan_direction)
        elif processing == 'differentiate':
            return self.differentiate(scan_direction)
        elif processing == 'subtract plane fit':
            return self.subtract_plane_fit(scan_direction)
        
    def raw (self, scan_direction):
        if scan_direction == 'fwd':
            z = self.signals['Z']['forward']
        elif scan_direction == 'bwd':
            import numpy as np
            z = np.flip (self.signals['Z']['backward'], axis = 1) # Reverse the order of elements in an array along the given axis. axis = 0: 행(상하) filp, axis = 1: 열(좌우) flip
        return z
    
    def subtract_average (self, scan_direction): # raw값에서 한 line의 평균값을 빼는것
        import warnings
        import numpy as np
        warnings.filterwarnings(action='ignore') # filterwarnings 무시
        z = self.raw(scan_direction)
        z_subav = np.zeros(np.shape(z)) # Return a new array of given shape and type, filled with zeros. z 와 같은크기의 array를 만듬
        lines = np.shape(z)[0] # z의 line수
        for i in range(lines):
            z_subav[i] = z[i] - np.nanmean(z[i]) # NaN을 무시하고 평균을 낸 값
        return z_subav

    def subtract_linear_fit (self, scan_direction): # raw값에서 linear fitting한 직선을 빼는것
        import numpy as np
        from scipy.optimize import curve_fit
        def f_lin(x, a, b): return a*x + b # curve fitting 형태 설정 (여기서는 linear)
        xrange = round(self.header['scan_range'][0] * 1e9)*1e-9 # x방향 scan range 조정 (0.xx nm 반올림) 
        z = self.raw(scan_direction)
        z_sublf = np.zeros(np.shape(z)) # z와 같은크기의 array를 만듬
        lines, pixels = np.shape(z) # z의 line 개수, line당 픽셀개수
        for i in range(lines):
            if np.shape(np.where(np.isnan(z))[0])[0] != 0: # image에 nan값이 포함되어 있을 경우 (== scan을 도중에 멈추었을 경우)
                if i < np.min(np.where(np.isnan(z))[0]): # scan한 부분까지
                    x = np.linspace(0, xrange, pixels)
                    popt, pcov = curve_fit(f_lin, x, z[i]) # x데이터, y데이터(=z[i])로 fitting 수행
                    z_sublf[i] = z[i] - f_lin(x, *popt)
                else:
                    z_sublf[i] = np.nan
            else:
                x = np.linspace(0, xrange, pixels)
                popt, pcov = curve_fit(f_lin, x, z[i]) # x - ith line: linear fitting
                z_sublf[i] = z[i] - f_lin(x, *popt)

        return z_sublf
    
    # 현재 이 함수는 line수와 pixel이 일치해야 사용가능
    def subtract_linear_fit_xy (self, scan_direction): # x,y축 모두에서 linear fitting 실행
        import numpy as np
        from scipy.optimize import curve_fit
        def f_lin(x, a, b): return a*x + b
        xrange = round(self.header['scan_range'][0] * 1e9)*1e-9
        z = self.subtract_linear_fit(scan_direction) # linear fit 실행
        z_sublf = np.zeros(np.shape(z))
        lines, pixels = np.shape(z)
        for i in range(lines):
            if np.shape(np.where(np.isnan(z))[0])[0] != 0: # image에 nan값이 포함되어 있을 경우 (== scan을 도중에 멈추었을 경우)
                if i < np.min(np.where(np.isnan(z))[0]):
                    x = np.linspace(0, xrange, pixels)
                    popt, pcov = curve_fit(f_lin, x, z.T[i]) # z를 transpose함 (== y축 방향으로 linear fit)
                    z_sublf[i] = z.T[i] - f_lin(x, *popt)
                else:
                    z_sublf[i] = np.nan
            else:
                x = np.linspace(0, xrange, pixels)
                popt, pcov = curve_fit(f_lin, x, z.T[i]) # x - ith line: linear fitting
                z_sublf[i] = z.T[i] - f_lin(x, *popt)

        return z_sublf.T # 다시 transpose 한 결과를 반환

    def subtract_parabolic_fit (self, scan_direction): # 2차함수로 fitting
        import numpy as np
        from scipy.optimize import curve_fit
        def f_parab(x, a, b, c): return a*(x**2) + b*x + c
        xrange = round(self.header['scan_range'][0] * 1e9)*1e-9
        z = self.raw(scan_direction)
        z_subpf = np.zeros(np.shape(z))
        lines, pixels = np.shape(z)
        for i in range(lines):
            if np.shape(np.where(np.isnan(z))[0])[0] != 0: # image에 nan값이 포함되어 있을 경우 (== scan을 도중에 멈추었을 경우)
                if i < np.min(np.where(np.isnan(z))[0]):
                    x = np.linspace(0, xrange, pixels)
                    popt, pcov = curve_fit(f_parab, x, z[i])
                    z_subpf[i] = z[i] - f_parab(x, *popt)
                else:
                    z_subpf[i] = np.nan
            else:
                x = np.linspace(0, xrange, pixels)
                popt, pcov = curve_fit(f_parab, x, z[i]) # x - ith line: linear fitting
                z_subpf[i] = z[i] - f_parab(x, *popt)

        return z_subpf
    
    def differentiate (self, scan_direction): # 미분값
        import numpy as np
        xrange, pixels = round(self.header['scan_range'][0] * 1e9)*1e-9, int(self.header['scan>pixels/line'])
        dx = xrange / pixels
        z = self.raw(scan_direction)
        z_deriv = np.zeros(np.shape(z))
        lines = np.shape(z)[0]
        for i in range(lines):
            z_deriv[i] = np.gradient(z[i], dx, edge_order = 2) # dI/dV curve를 직접 미분. --> d^2I/dV^2
        return z_deriv
    
    def subtract_plane_fit (self, scan_direction):
        import numpy as np
        from scipy.linalg import lstsq
        
        # regular grid covering the domain of the data
        Z = self.raw(scan_direction)
        X, Y = np.meshgrid( np.arange(np.shape(Z)[1]), np.arange(np.shape(Z)[0]) ) # x-y plain에 point 균등하게 찍기 (plane 생성)
        
        # best-fit linear plane
        A = np.c_[X.flatten(), Y.flatten(), np.ones( np.shape(Z.flatten())[0] )]
        C, _, _, _ = lstsq(A, Z.flatten())    # coefficients

        # evaluate it on grid
        plane = C[0]*X + C[1]*Y + C[2]
        return Z - plane

class didvmap:
    
    '''
    Args:
        filepath : str
            Name of the Nanonis spectrum file to be loaded.
    
    Attributes (name : type):
        fname : str
            The name of the file excluding its containing directory.
        header : dict
            Header information of spectrum data.
        signals : dict
            Measured values in spectrum data.
    
    Methods:
        get_map(self, scan_direction = 'fwd', channel = 'LI_Demod_1_X')
            Parameters:
            scan_direction : str
                The direction of scan.
                Possible parameters is the following: 
                    'fwd', 'bwd'
            channel : str
                The channel to be returned.
                'LI_Demod_1_X' channel is returned by default.
                Other channels can also be returned if the input file contains. e.g. 'LI_Demod_1_Y', 'LI_Demod_2_X', ...
            
            Returns the two-dimensional array
            that represents the dI/dV map data scanned in the direction of scan_direction.
    '''
    
    def __init__(self, filepath):
        import nanonispy as nap
        self.fname = nap.read.Scan(filepath).fname.split('\\')[-1]
        self.header = nap.read.Scan(filepath).header
        self.signals = nap.read.Scan(filepath).signals
    
    def get_map(self, scan_direction = 'fwd', channel = 'LI_Demod_1_X'): # dI/dV 채널에서 가져오기
        if scan_direction == 'fwd':
            didv = self.signals[channel]['forward']
        elif scan_direction == 'bwd':
            import numpy as np
            didv = np.flip(self.signals[channel]['backward'], axis = 1)
        return didv

class currentmap:
    
    '''
    Args:
        filepath : str
            Name of the Nanonis spectrum file to be loaded.
    
    Attributes (name : type):
        fname : str
            The name of the file excluding its containing directory.
        header : dict
            Header information of spectrum data.
        signals : dict
            Measured values in spectrum data.
    
    Methods:
        get_map(self, scan_direction = 'fwd')
            Parameters:
            scan_direction : str
                The direction of scan.
                Possible parameters is the following: 
                    'fwd', 'bwd'
            
            Returns the two-dimensional array
            that represents the current map data scanned in the direction of scan_direction.
    '''
    
    def __init__(self, filepath):
        import nanonispy as nap
        self.fname = nap.read.Scan(filepath).fname.split('\\')[-1]
        self.header = nap.read.Scan(filepath).header
        self.signals = nap.read.Scan(filepath).signals
    
    def get_map(self, scan_direction = 'fwd'):
        if scan_direction == 'fwd':
            current = self.signals['Current']['forward']
        elif scan_direction == 'bwd':
            import numpy as np
            current = np.flip(self.signals['Current']['backward'], axis = 1)
        return current

class fft:
    
    '''
    Args:
        filepath : str
            Name of the Nanonis spectrum file to be loaded.
    
    Attributes (name : type):
        fname : str
            The name of the file excluding its containing directory.
        header : dict
            Header information of spectrum data.
        signals : dict
            Measured values in spectrum data.
    '''
    def __init__(self, filepath):
        import nanonispy as nap
        self.fname = nap.read.Scan(filepath).fname.split('/')[-1]
        self.header = nap.read.Scan(filepath).header
        self.signals = nap.read.Scan(filepath).signals
        
    def two_d_FFT_sqrt(image):
        '''
        Parameters
        ----------
        image : 2D numpy input
        Calculate FFT

        Returns
        -------
        image_fit : 2D numpy output
        Check the scale bar
        '''
        import numpy as np
        fft = np.fft.fft2(image) # FFT only
        fft_shift = np.fft.fftshift(fft) # 2D FFT 를 위한 image shift, Shift the zero-frequency component to the center of the spectrum.
        image_fft = np.sqrt(np.abs(fft_shift)) # 절댓값의 제곱근
        return image_fft
    
    def two_d_FFT_log(image):
        '''
        Parameters
        ----------
        image : 2D numpy input
        Calculate FFT

        Returns
        -------
        image_fit : 2D numpy output
        Check the scale bar
        '''
        import numpy as np
        fft = np.fft.fft2(image) # FFT only
        fft_shift = np.fft.fftshift(fft) #2D FFT 를 위한 image shift
        image_fft = np.log(np.abs(fft_shift)) # 절댓값의 로그
        return image_fft
    
    def two_d_FFT_lin(image):
        '''
        Parameters
        ----------
        image : 2D numpy input
        Calculate FFT

        Returns
        -------
        image_fit : 2D numpy output
        Check the scale bar
        '''
        import numpy as np
        fft = np.fft.fft2(image) # FFT only
        fft_shift = np.fft.fftshift(fft) #2D FFT 를 위한 image shift
        image_fft = np.abs(fft_shift) # 절댓값
        return image_fft
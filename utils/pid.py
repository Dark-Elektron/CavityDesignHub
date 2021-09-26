from termcolor import colored


class PID:
    def __init__(self, kp, ki, kd, setpoint, limits=None):
        if limits is None:
            self.limits = [0, 1]
        else:
            self.limits = limits

        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.setpoint = setpoint
        self.error = 0
        self.integral_error = 0
        self.error_last = 0
        self.derivative_error = 0
        self.output = 0
        self.dt = 0
        self.count = 0

    def compute(self, currentpoint, dt):
        self.currentpoint = currentpoint
        self.error = self.setpoint - currentpoint
        self.integral_error += self.error * dt
        self.derivative_error = (self.error - self.error_last)/dt

        self.output = self.kp*self.error + self.ki*self.integral_error + self.kd*self.derivative_error

        self.error_last = self.error
        print(colored(f'\t\tOutput, error: {self.output, self.error}', 'yellow'))
        print()
        if self.output >= self.limits[1]:
            return self.limits[1] #, self.error
        elif self.output <= self.limits[0]:
            return self.limits[0] #, self.error
        else:
            return self.output #, self.error

    def setSetpoint(self, sp):
        self.setpoint = sp

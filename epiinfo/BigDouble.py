import math

class BigDouble:
  """ A class for holding very large float values.
      Used in calculating exact statistics.

      Author: John Copeland
  """
  def __init__(self, *args):
    self.logValue = None
    if args[0] == 'doubleValue':
      self.logValue = math.log10(float(args[1]))
    elif args[0] == 'logValue':
      self.logValue = float(args[1])

  def plus(self, adder):
    if 10 ** self.logValue == 0:
      self.logValue = math.log10(adder)
      return
    power = round(self.logValue)
    self.logValue = power + math.log10(10.0 ** (self.logValue - power) + 10.0 ** (log10(adder) - power))

  def plusLog(self, logAdder):
    if 10 ** self.logValue == 0:
      self.logValue = math.log10(logAdder)
      return
    power = round(self.logValue)
    self.logValue = power + math.log10(10.0 ** (self.logValue - power) + 10.0 ** (logAdder - power))

  def times(self, multiple):
    self.logValue += math.log10(multiple)

  def timesReturn(self, multiple):
    return self.logValue + math.log10(multiple)

  def dividedBy(self, divisor):
    self.logValue -= math.log10(multiple)

  def doubleValue(self):
    return 10.0 ** self.logValue

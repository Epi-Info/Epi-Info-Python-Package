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
    if adder > 0:
      self.logValue = power + math.log10(10.0 ** (self.logValue - power) + 10.0 ** (math.log10(adder) - power))
    elif adder < 0:
      try:
        self.logValue = power + math.log10(10.0 ** (self.logValue - power) - 10.0 ** (math.log10(-adder) - power))
      except ValueError:
        self.logValue = 0.0

  def plusLog(self, logAdder):
    if 10 ** self.logValue == 0:
      self.logValue = logAdder
      return
    power = round(self.logValue)
    self.logValue = power + math.log10(10.0 ** (self.logValue - power) + 10.0 ** (logAdder - power))

  def times(self, multiple):
    try:
      self.logValue += math.log10(multiple)
    except ValueError:
      self.logValue = float('-inf')

  def timesReturn(self, multiple):
    try:
      return self.logValue + math.log10(multiple)
    except ValueError:
      return float('-inf')

  def dividedBy(self, divisor):
    try:
      self.logValue -= math.log10(divisor)
    except ValueError:
      self.logValue = float('inf')

  def doubleValue(self):
    return 10.0 ** self.logValue

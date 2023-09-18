class Series:
	def __init__(self, start_year, end_year):
		self._list = []
		self.start_year = start_year
		self.end_year = end_year
		
		for a in range(end_year - start_year + 1):
			self._list.append(None)
	
	def __len__(self):
		return len(self._list)
	
	def __getitem__(self, year):
		if isinstance(year, slice):
			if year.start < self.start_year:
				raise IndexError
			return self._list[year.start - self.start_year:year.stop - self.start_year + 1]
		else:
			if year < self.start_year or year > self.end_year:
				raise IndexError
			return self._list[year - self.start_year]
	
	def __setitem__(self, year, value):
		if isinstance(year, slice):
			if year.start < self.start_year:
				raise IndexError
			if not isinstance(value, list):
				value = [value]*(year.stop-year.start+1)
			self._list[year.start - self.start_year:year.stop - self.start_year + 1] = value
		else:
			if year < self.start_year or year > self.end_year:			
				raise IndexError
			self._list[year - self.start_year] = value
	
